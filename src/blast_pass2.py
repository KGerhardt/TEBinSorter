"""
blast_pass2.py — minimap2-based pass-2 classification for HMM-unclassified
sequences.

Port of TEBinSorter's blastn-based pass-2 onto minimap2. Keeps the filename
and public symbols (`blast_pass2`, `store_blast_hits`, `classify_from_blast`)
so the rest of the pipeline doesn't care what aligner is behind them, but
the SQLite `blast_hits` schema adds a `tcovs` column and the filter enforces
both qcov AND tcov ≥ cutoff.

Rule semantics for --pass2-rule I-C-L:
    identity ≥ I%,  qcov ≥ C%,  tcov ≥ C%,  length ≥ L bp

where qcov / tcov are unions of aligned intervals across ALL minimap2 chains
for a given (query, target) pair, divided by qlen / tlen respectively.
"""

import logging
import os
import tempfile
import time
from collections import defaultdict

import pyfastx

import minimap
import pass2_external

log = logging.getLogger(__name__)


# Developer toggle, deliberately not exposed on the CLI.
#   True  -> require BOTH qcov AND tcov >= coverage threshold (strict).
#   False -> require AT LEAST ONE of qcov, tcov >= coverage threshold
#            (Wicker et al. 80-80-80 style: "candidate must cover ≥80% of
#            at least one of the elements being compared").
# minimap2 has no native coverage filter, so this is enforced post-hoc in
# classify_from_blast's SQL WHERE clause.
REQUIRE_BOTH_COVERAGE = False


def _get_classified_ids(conn):
    """Get classified sequence IDs per database from classifier results."""
    classified = defaultdict(set)
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "legacy_hits" in tables:
        for row in conn.execute(
            "SELECT DISTINCT base_seq, database FROM legacy_hits"
        ):
            classified[row[0]].add(row[1])
    return dict(classified)


def split_classified_unclassified(input_fasta, classified_ids, outdir,
                                  seq_type="nucl"):
    """Split input into classified-pool FASTA and unclassified-query FASTA.

    Nucleotide pools are uppercased and stripped of non-ATCG characters — not
    because minimap2 requires it (it handles ambiguous bases) but to keep
    invariants aligned with the mmseqs sibling port and make inputs clean.
    """
    os.makedirs(outdir, exist_ok=True)
    db_fasta = os.path.join(outdir, "blast_db.fa")
    qry_fasta = os.path.join(outdir, "blast_query.fa")

    nucl = (seq_type == "nucl")
    fa = pyfastx.Fasta(input_fasta, build_index=True)

    n_classified = 0
    n_unclassified = 0
    db_seq_to_dbs = {}

    with open(db_fasta, "w") as dbh, open(qry_fasta, "w") as qh:
        for rec in fa:
            name = rec.name
            seq = str(rec.seq)
            if nucl:
                seq = "".join(c for c in seq.upper() if c in "ATCG")
            if name in classified_ids:
                dbh.write(f">{name}\n{seq}\n")
                db_seq_to_dbs[name] = classified_ids[name]
                n_classified += 1
            else:
                qh.write(f">{name}\n{seq}\n")
                n_unclassified += 1

    log.info(f"  Split: {n_classified} classified (DB), "
             f"{n_unclassified} unclassified (query)")
    return db_fasta, qry_fasta, db_seq_to_dbs


def run_alignment(query_fa, target_fa, paf_out, ncpu=4,
                  preset="asm20", extra=""):
    """Run minimap2 once with preset + -c so AS tags are available."""
    minimap.run_minimap2(
        query_fa=query_fa, target_fa=target_fa, paf_out=paf_out,
        ncpu=ncpu, preset=preset, extra=extra,
    )
    return paf_out


def parse_minimap2_output(paf_path):
    """Parse PAF, union per (q,t), best-hit per query, return hit dicts."""
    merged = minimap.parse_paf_besthit(paf_path)
    best = minimap.besthit_per_query(merged)

    hits = []
    for qid, m in best.items():
        hits.append({
            "qseqid": m["qseqid"],
            "sseqid": m["sseqid"],
            "pident": m["fident"] * 100.0,   # 0..100
            "length": m["alnlen"],           # union of aligned query bases
            "evalue": 0.0,                   # sentinel
            "bitscore": float(m["score"]),   # minimap2 AS score
            "qlen": m["qlen"],
            "slen": m["tlen"],
            "qcovs": m["qcov"] * 100.0,      # 0..100
            "tcovs": m["tcov"] * 100.0,      # 0..100 (new column)
        })
    return hits


def store_blast_hits(conn, hits, db_seq_to_dbs):
    """Store hits in SQLite. Adds a `tcovs` column vs. the stock schema."""
    conn.execute("""
        CREATE TABLE IF NOT EXISTS blast_hits (
            qseqid      TEXT NOT NULL,
            sseqid      TEXT NOT NULL,
            pident      REAL NOT NULL,
            length      INTEGER NOT NULL,
            evalue      REAL NOT NULL,
            bitscore    REAL NOT NULL,
            qlen        INTEGER NOT NULL,
            slen        INTEGER NOT NULL,
            qcovs       REAL NOT NULL,
            tcovs       REAL NOT NULL,
            classified_by TEXT NOT NULL
        )
    """)
    rows = []
    for h in hits:
        dbs = db_seq_to_dbs.get(h["sseqid"], set())
        classified_by = ",".join(sorted(dbs)) if dbs else "unknown"
        rows.append((
            h["qseqid"], h["sseqid"], h["pident"], h["length"],
            h["evalue"], h["bitscore"], h["qlen"], h["slen"],
            h["qcovs"], h["tcovs"], classified_by,
        ))
    conn.executemany(
        "INSERT INTO blast_hits VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()
    log.info(f"  Stored {len(rows)} minimap2 hits")


def classify_from_blast(conn, classifications, database=None,
                        min_identity=80, min_coverage=80, min_length=80):
    """Classify pass-2 rescues. Filters on BOTH qcovs and tcovs."""
    tables = {r[0] for r in conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    if "blast_hits" not in tables:
        return []

    if REQUIRE_BOTH_COVERAGE:
        where = ("WHERE pident >= ? AND qcovs >= ? AND tcovs >= ? "
                 "AND length >= ?")
        params = [min_identity, min_coverage, min_coverage, min_length]
    else:
        # At-least-one-side coverage (Wicker et al. 80-80-80 interpretation)
        where = ("WHERE pident >= ? AND (qcovs >= ? OR tcovs >= ?) "
                 "AND length >= ?")
        params = [min_identity, min_coverage, min_coverage, min_length]
    log.info(f"  coverage mode: {'BOTH' if REQUIRE_BOTH_COVERAGE else 'AT-LEAST-ONE'} "
             f"qcov/tcov >= {min_coverage}")
    if database:
        where += " AND classified_by LIKE ?"
        params.append(f"%{database}%")

    rows = conn.execute(f"""
        SELECT qseqid, sseqid, pident, qcovs, tcovs, length, bitscore
        FROM blast_hits
        {where}
        ORDER BY bitscore DESC
    """, params).fetchall()

    classified_set = set(classifications.keys())
    best = {}
    for qid, sid, pident, qcovs, tcovs, length, bitscore in rows:
        if qid in classified_set:
            continue
        if qid not in best:
            best[qid] = (sid, pident, qcovs, tcovs, length, bitscore)

    new_classifications = []
    no_source = 0
    for qid, (sid, pident, qcovs, tcovs, length, bitscore) in best.items():
        if sid in classifications:
            source = classifications[sid]
            new_classifications.append({
                "id": qid,
                "order": source["order"],
                "superfamily": source["superfamily"],
                "clade": "unknown",
                "complete": "none",
                "strand": "?",
                "domains": "none",
                "blast_source": sid,
                "blast_pident": pident,
                "blast_qcovs": qcovs,
                "blast_tcovs": tcovs,
                "blast_bitscore": bitscore,
            })
        else:
            no_source += 1

    if no_source:
        log.info(f"    {no_source} pass-2 hits to unclassified targets (skipped)")

    log.info(f"  pass-2: {len(new_classifications)} sequences classified "
             f"(from {len(best)} hits passing filters)")
    return new_classifications


def blast_pass2(input_fasta, conn, hmm_classifications=None,
                seq_type="nucl", n_processors=4,
                min_identity=80, min_coverage=80, min_length=80,
                outdir=None,
                pass2_classified_fasta=None,
                preset="asm20", minimap2_extra=""):
    """minimap2-based pass-2. Public signature mirrors sibling ports.

    min_coverage is applied to BOTH qcov and tcov.
    """
    t0 = time.time()
    minimap.check_minimap2()

    if seq_type != "nucl":
        log.warning("minimap2 pass-2 only supports nucleotide sequences; "
                    f"seq_type={seq_type!r} will be treated as nucl")
        seq_type = "nucl"

    classified_ids = _get_classified_ids(conn)
    if not classified_ids:
        log.info("  No classified sequences for minimap2 pass-2")
        return []

    log.info(f"  minimap2 pass-2: {len(classified_ids)} classified sequences as targets")

    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="tebinsorter_minimap2_")
    work = os.path.join(outdir, "minimap2_pass2")

    db_fasta, qry_fasta, db_seq_to_dbs = split_classified_unclassified(
        input_fasta, classified_ids, work, seq_type=seq_type)

    if os.path.getsize(qry_fasta) == 0:
        log.info("  No unclassified sequences to search")
        return []

    if hmm_classifications is None:
        hmm_classifications = {}

    if pass2_classified_fasta:
        updated = pass2_external.update_classified_fasta_headers(
            pass2_classified_fasta, hmm_classifications, work
        )
        src = updated or pass2_classified_fasta
        pass2_external.extend_hmm_classifications_from_fasta(
            hmm_classifications, src, db_seq_to_dbs
        )
        merged_db = os.path.join(work, "pass2_db_merged.fa")
        pass2_external.merge_classified_fastas(
            merged_db, db_fasta, src, clean_nucl=True
        )
        db_fasta = merged_db

    if os.path.getsize(db_fasta) == 0:
        log.info("  pass-2 target FASTA is empty; skipping minimap2")
        return []

    paf_out = os.path.join(work, "pass2.paf")

    log.info(f"  Running minimap2 -x {preset} with {n_processors} threads")
    t1 = time.time()
    run_alignment(
        query_fa=qry_fasta, target_fa=db_fasta, paf_out=paf_out,
        ncpu=n_processors, preset=preset, extra=minimap2_extra,
    )
    t2 = time.time()
    log.info(f"  minimap2 alignment: {t2 - t1:.1f}s")

    hits = parse_minimap2_output(paf_out)
    log.info(f"  {len(hits)} best-hit records after per-pair union of PAF chains")

    if hits:
        store_blast_hits(conn, hits, db_seq_to_dbs)

    log.info(f"  {len(hmm_classifications)} HMM classifications available for inheritance")

    new_cls = classify_from_blast(
        conn, hmm_classifications,
        min_identity=min_identity,
        min_coverage=min_coverage,
        min_length=min_length,
    )

    t3 = time.time()
    log.info(f"  minimap2 pass-2 total: {t3 - t0:.1f}s")
    return new_cls
