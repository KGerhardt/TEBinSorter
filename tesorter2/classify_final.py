#!/usr/bin/env python3
"""Classify minimap2 PAF queries at mutual-coverage + identity thresholds.

Auto-adapts to PAF flavor:
  * with -c   : has cg:Z, de:f, NM:i tags + exact col10/col11. Uses de:f.
  * without -c: only dv:f (chain-estimated divergence) + approximate
                col10 (= k * num_seeds) and col11 (= span). Uses dv:f.
  * legacy    : neither tag — falls back to matches/block_len (BLAST-style),
                which over-penalises structural indels (a single 3 kb
                deletion drops identity from ~99% to ~63%).

Per (query,target) pair:
  - drops alignments whose query span is fully encompassed by a strictly
    longer one (avoids double-counting overlapped regions)
  - identity = qspan-weighted mean of (1 - div) over surviving rows
  - qcov / tcov = union-of-intervals length divided by qlen / tlen

Output (TSV):
  query  qlen  best_target  tlen  qcov  tcov  identity  ident_src  pass
ident_src in {de, dv, mb}.
"""
import argparse
import sys
from collections import defaultdict


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for a, b in ivs[1:]:
        if a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def total_len(ivs):
    return sum(b - a for a, b in ivs)


def parse_row(line):
    """Return dict for an aligned row, ('NOHIT', qname, qlen) for paf-no-hit
    rows (target == '*'), or None for malformed rows."""
    f = line.rstrip("\n").split("\t")
    if len(f) < 12:
        return None
    qname = f[0]
    try:
        qlen = int(f[1])
    except ValueError:
        return None
    if f[5] == "*":
        return ("NOHIT", qname, qlen)
    try:
        rec = {
            "q":      qname,
            "qlen":   qlen,
            "qs":     int(f[2]),
            "qe":     int(f[3]),
            "t":      f[5],
            "tlen":   int(f[6]),
            "ts":     int(f[7]),
            "te":     int(f[8]),
            "matches": int(f[9]),
            "blocks":  int(f[10]),
        }
    except ValueError:
        return None
    de = dv = None
    for tag in f[12:]:
        if tag.startswith("de:f:"):
            try:
                de = float(tag[5:])
            except ValueError:
                pass
        elif tag.startswith("dv:f:"):
            try:
                dv = float(tag[5:])
            except ValueError:
                pass
    rec["de"] = de
    rec["dv"] = dv
    return rec


def aggregate_pair(rows, qlen, tlen):
    """Return (qcov, tcov, identity, ident_src) for one (q,t) pair."""
    # Coverage from union of spans
    qcov = (total_len(merge_intervals([(r["qs"], r["qe"]) for r in rows]))
            / qlen) if qlen else 0.0
    tcov = (total_len(merge_intervals([(r["ts"], r["te"]) for r in rows]))
            / tlen) if tlen else 0.0

    # Encompass-drop on query: discard alignments whose qspan is fully
    # contained in a strictly longer row's qspan
    survivors = []
    for r in rows:
        rspan = r["qe"] - r["qs"]
        encompassed = False
        for o in rows:
            if o is r:
                continue
            ospan = o["qe"] - o["qs"]
            if ospan > rspan and o["qs"] <= r["qs"] and r["qe"] <= o["qe"]:
                encompassed = True
                break
        if not encompassed:
            survivors.append(r)
    if not survivors:
        survivors = rows  # degenerate, but keep something

    # Identity source priority: de:f > dv:f > matches/blocks
    has_de = any(r["de"] is not None for r in survivors)
    has_dv = any(r["dv"] is not None for r in survivors)

    if has_de:
        wsum = 0.0; wtot = 0
        for r in survivors:
            if r["de"] is None:
                continue
            w = r["qe"] - r["qs"]
            wsum += r["de"] * w
            wtot += w
        ident = 1.0 - (wsum / wtot if wtot else 1.0)
        src = "de"
    elif has_dv:
        wsum = 0.0; wtot = 0
        for r in survivors:
            if r["dv"] is None:
                continue
            w = r["qe"] - r["qs"]
            wsum += r["dv"] * w
            wtot += w
        ident = 1.0 - (wsum / wtot if wtot else 1.0)
        src = "dv"
    else:
        m = sum(r["matches"] for r in rows)
        b = sum(r["blocks"] for r in rows)
        ident = (m / b) if b else 0.0
        src = "mb"

    return qcov, tcov, ident, src


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paf", help="input PAF (use '-' for stdin)")
    ap.add_argument("-o", "--out", default="-",
                    help="output TSV (default stdout)")
    ap.add_argument("--qcov",  type=float, default=0.70)
    ap.add_argument("--tcov",  type=float, default=0.70)
    ap.add_argument("--ident", type=float, default=0.70)
    ap.add_argument("--queries-fa", default=None,
                    help="optional FASTA so unmapped queries also appear")
    args = ap.parse_args()

    pairs = defaultdict(list)
    pair_lens = {}
    qlen_seen = {}
    fmt_de = fmt_dv = fmt_mb = 0

    fh = sys.stdin if args.paf == "-" else open(args.paf)
    for line in fh:
        rec = parse_row(line)
        if rec is None:
            continue
        if isinstance(rec, tuple) and rec[0] == "NOHIT":
            qlen_seen[rec[1]] = rec[2]
            continue
        qlen_seen[rec["q"]] = rec["qlen"]
        key = (rec["q"], rec["t"])
        pairs[key].append(rec)
        pair_lens[key] = (rec["qlen"], rec["tlen"])
        if rec["de"] is not None:
            fmt_de += 1
        elif rec["dv"] is not None:
            fmt_dv += 1
        else:
            fmt_mb += 1

    sys.stderr.write(
        f"[classify_final] parsed pairs: {len(pairs)}; "
        f"rows w/ de:f={fmt_de}, dv:f={fmt_dv}, neither={fmt_mb}\n"
    )

    # Best (q,t) per query
    best = {}
    for (q, t), rows in pairs.items():
        qlen, tlen = pair_lens[(q, t)]
        qcov, tcov, ident, src = aggregate_pair(rows, qlen, tlen)
        passed = qcov >= args.qcov and tcov >= args.tcov and ident >= args.ident
        score = (1 if passed else 0, ident * min(qcov, tcov))
        cur = best.get(q)
        if cur is None or score > cur[0]:
            best[q] = (score, qlen, t, tlen, qcov, tcov, ident, src, passed)

    # Optional FASTA enumeration to fill in unmapped queries
    if args.queries_fa:
        with open(args.queries_fa) as fa:
            cur_id = None; cur_len = 0
            for ln in fa:
                if ln.startswith(">"):
                    if cur_id is not None and cur_id not in best:
                        best[cur_id] = ((0, 0), cur_len, "*", 0,
                                        0.0, 0.0, 0.0, "-", False)
                    cur_id = ln[1:].split()[0]
                    cur_len = 0
                else:
                    cur_len += len(ln.strip())
            if cur_id is not None and cur_id not in best:
                best[cur_id] = ((0, 0), cur_len, "*", 0,
                                0.0, 0.0, 0.0, "-", False)
    else:
        # at least include queries seen in PAF (mapped + paf-no-hit)
        for q, qlen in qlen_seen.items():
            if q not in best:
                best[q] = ((0, 0), qlen, "*", 0, 0.0, 0.0, 0.0, "-", False)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write("query\tqlen\tbest_target\ttlen\tqcov\ttcov\t"
              "identity\tident_src\tpass\n")
    n_pass = 0
    for q in sorted(best):
        (_, qlen, t, tlen, qcov, tcov, ident, src, passed) = best[q]
        if passed:
            n_pass += 1
        out.write(f"{q}\t{qlen}\t{t}\t{tlen}\t{qcov:.4f}\t{tcov:.4f}"
                  f"\t{ident:.4f}\t{src}\t{int(passed)}\n")

    sys.stderr.write(
        f"[classify_final] queries: {len(best)}  "
        f"classified: {n_pass}\n"
    )


if __name__ == "__main__":
    main()
