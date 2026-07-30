"""Assert that a multi-database TEsorter2 run harmonizes clade nomenclature.

Reads the combined ``{prefix}.cls.tsv`` from a ``-d rexdb,gydb`` run and checks
that the cross-database vote ran on harmonized labels rather than raw strings.

The substantive check is the concordance contrast: REXdb and GyDB share almost
no clade spellings, so agreement measured on raw tokens is near zero even where
the two databases name the same lineage. Harmonized agreement must be strictly
higher, and the combined call must never be left on a GyDB-only native token.

Usage: check_multidb.py <combined.cls.tsv>
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from tesorter2.clade_harmonize import (load_harmonization, harmonize,  # noqa: E402
                                       harmonize_lineage)

fails = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))
    if not ok:
        fails.append(name)


def parse_secondary(field):
    """{db: (order, superfamily, clade, score)} from the SecondaryHits column.

    Native GyDB clades can themselves contain '/' (e.g. ``ty1/copia``), so the
    triple is split on the first two separators only.
    """
    out = {}
    if not field or field == ".":
        return out
    for hit in field.split(";"):
        db, _, rest = hit.partition(":")
        triple, _, score = rest.rpartition("=")
        parts = triple.split("/", 2)
        if len(parts) != 3:
            continue
        out[db] = (parts[0], parts[1], parts[2], float(score))
    return out


def main(path):
    H = load_harmonization()
    check("harmonization table is loaded", len(H) > 0, f"{len(H)} rows")

    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
    check("SecondaryHits column present", "SecondaryHits" in header)
    if "SecondaryHits" not in header:
        return 1
    sec_i = header.index("SecondaryHits")

    # GyDB native clade tokens, and the harmonized names no database emits
    # natively. Either set appearing on the wrong side of the run is decisive.
    gydb_native = {c for (db, _o, _s, c) in H if db == "gydb"}
    unified_only = {v[2] for k, v in H.items() if k[0] == "gydb"} - gydb_native

    n_multi = raw_agree = harm_agree = 0
    leaked = []
    off_taxonomy = []
    unified_seen = set()
    native_secondary = 0

    for r in rows:
        clade = r[3]
        sec = parse_secondary(r[sec_i])
        if clade in unified_only:
            unified_seen.add(clade)
        if "gydb" not in sec or "rexdb" not in sec:
            continue
        n_multi += 1

        g_o, g_sf, g_c, _ = sec["gydb"]
        x_o, x_sf, x_c, _ = sec["rexdb"]
        if g_c in gydb_native:
            native_secondary += 1
        if g_c == x_c:
            raw_agree += 1
        if harmonize("gydb", g_o, g_sf, g_c, H)[2] == \
           harmonize("rexdb", x_o, x_sf, x_c, H)[2]:
            harm_agree += 1
        # The combined call must not be a GyDB-only lowercase native token:
        # that is exactly what an un-harmonized vote leaves behind.
        if clade in gydb_native and clade not in unified_only:
            leaked.append((r[0], clade))
        # Stronger, and symmetric across databases: the reported clade must be
        # one of the *harmonized* candidate labels. An un-harmonized vote
        # reports the winner's native token, which generally is not one of them
        # (REXdb 'Ale' vs the harmonized 'Ale/Retrofit'), so this catches a
        # regression on either side rather than only on GyDB's.
        candidates = {harmonize(db, o, sf, c, H)[2]
                      for db, (o, sf, c, _s) in sec.items()}
        if clade not in candidates:
            off_taxonomy.append((r[0], clade, sorted(candidates)))

    print(f"\n  elements: {len(rows)}   classified by both databases: {n_multi}")
    if n_multi:
        print(f"  clade agreement, raw tokens : {raw_agree}/{n_multi} "
              f"({100.0 * raw_agree / n_multi:.1f}%)")
        print(f"  clade agreement, harmonized : {harm_agree}/{n_multi} "
              f"({100.0 * harm_agree / n_multi:.1f}%)\n")

    # The clade name must not depend on which database won the score vote.
    # Every element whose databases fall in the same group must report the same
    # Clade, with the finer call (if any) in Lineage.
    if "Lineage" in header:
        lin_i = header.index("Lineage")
        by_pair, bad = {}, []
        for r in rows:
            sec = parse_secondary(r[sec_i])
            if "gydb" not in sec or "rexdb" not in sec:
                continue
            key = (sec["gydb"][2], sec["rexdb"][2])   # native pair
            by_pair.setdefault(key, set()).add(r[3])
            # A database that resolves a sub-lineage must not have it dropped.
            want = {harmonize_lineage(db, o, sf, c, H)
                    for db, (o, sf, c, _s) in sec.items()} - {""}
            if want and r[lin_i] == ".":
                bad.append((r[0], sorted(want)))
        forked = {k: sorted(v) for k, v in by_pair.items() if len(v) > 1}
        check("same native pair always yields the same Clade", not forked,
              f"{len(forked)} forked: {list(forked.items())[:3]}")
        check("a resolved sub-lineage is never dropped", not bad,
              f"{len(bad)} dropped: {bad[:3]}")

    check("both databases classified a shared set of elements", n_multi > 0,
          f"n={n_multi}")
    check("SecondaryHits preserves native GyDB spellings", native_secondary > 0,
          f"{native_secondary} rows")
    check("harmonization raises cross-database clade agreement",
          harm_agree > raw_agree, f"{raw_agree} -> {harm_agree}")
    check("combined call never left on a GyDB-only native token",
          not leaked, f"{len(leaked)} leaked: {leaked[:5]}")
    check("combined clade is always a harmonized candidate label",
          not off_taxonomy, f"{len(off_taxonomy)} off: {off_taxonomy[:3]}")
    if unified_seen:
        print(f"  unified names in the combined call: {sorted(unified_seen)}")

    print()
    if fails:
        print(f"RESULT: FAIL - {len(fails)} check(s) failed: {fails}")
        return 1
    print("RESULT: PASS - multi-database run is harmonized")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
