"""
Emit routed FASTA partitions for external aligners (e.g. BATH).

Takes pass-1 results from the database, plans the same work units
pass-2 would have used, and writes labeled FASTA files:
    {MODEL_NAME}_partition_{BATCH}.fasta

Each file contains the sequences routed to that model in that batch.
Sequences may appear in multiple files if they route to multiple models.

Can be called standalone:
    python3 emit.py results.db sequences.faa hmm_db.hmm -o output.BATHwater
"""

import os
import sys
import logging
import argparse
import sqlite3
from collections import defaultdict

from sequence import load_sequences_dict, parse_frame_suffix
from search import _invert_seq_models, _plan_work_units
from hmm import load_hmms

log = logging.getLogger(__name__)


def _sanitize_model_name(model_name):
    """Make a model name safe for use as a filename."""
    return model_name.replace("/", "_").replace(":", "_").replace(" ", "_")


def emit_partitions(conn, seq_fasta, hmm_path, output_dir, n_workers=4,
                    db_name=None):
    """
    Emit routed FASTA partitions from pass-1 results.

    Args:
        conn: sqlite3 connection with pass1_hits populated
        seq_fasta: path to the sequence FASTA (translated AA or raw nucl)
        hmm_path: path to HMM database file (for M values)
        output_dir: directory to write partition FASTA files
        n_workers: number of workers (controls partition granularity)
        db_name: if set, filter to this database only

    Returns:
        list of (model_name, fasta_path) tuples
    """
    os.makedirs(output_dir, exist_ok=True)

    # Rebuild seq_models from pass1_hits
    query = "SELECT target_name, query_name FROM pass1_hits"
    params = ()
    if db_name is not None:
        query += " WHERE database = ?"
        params = (db_name,)

    seq_models = defaultdict(set)
    for row in conn.execute(query, params):
        seq_models[row[0]].add(row[1])
    seq_models = dict(seq_models)

    if not seq_models:
        log.info("No pass-1 hits to emit")
        return []

    # Load HMMs for M values
    hmms = load_hmms(hmm_path)
    hmms_dict = {h.name: h for h in hmms}

    # Load sequence lengths
    import pyfastx
    fa = pyfastx.Fasta(seq_fasta, build_index=True)
    seq_lens = {rec.name: len(rec) for rec in fa}

    # Invert and plan work units
    model_seqs = _invert_seq_models(seq_models)
    model_seqs = {m: s for m, s in model_seqs.items() if m in hmms_dict}
    work_units = _plan_work_units(model_seqs, hmms_dict, seq_lens, n_workers)

    if not work_units:
        log.info("No work units to emit")
        return []

    # Load all sequences into memory for fast slicing
    seqs = load_sequences_dict(seq_fasta)

    # Track batch numbers per model
    model_batch_counts = defaultdict(int)
    emitted = []

    for model_names, seq_names, cost in work_units:
        for model_name in model_names:
            batch_num = model_batch_counts[model_name]
            model_batch_counts[model_name] += 1

            safe_name = _sanitize_model_name(model_name)
            fname = f"{safe_name}_partition_{batch_num}.fasta"
            fpath = os.path.join(output_dir, fname)

            # Only emit sequences that route to this specific model
            model_seq_set = model_seqs.get(model_name, set())
            batch_seqs = set(seq_names) & model_seq_set

            if not batch_seqs:
                continue

            with open(fpath, "w") as f:
                for seq_name in sorted(batch_seqs):
                    if seq_name in seqs:
                        f.write(f">{seq_name}\n{seqs[seq_name]}\n")

            emitted.append((model_name, fpath))
            log.info(f"  {fname}: {len(batch_seqs)} sequences")

    log.info(f"Emitted {len(emitted)} partition files to {output_dir}")
    return emitted


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="emit-bath",
        description="Emit routed FASTA partitions from a TEBinSorter "
                    "results database for use with the BATH aligner.",
    )
    parser.add_argument(
        "database",
        help="Path to TEBinSorter results .db file",
    )
    parser.add_argument(
        "sequences",
        help="Path to the sequence FASTA (translated AA or raw nucleotide)",
    )
    parser.add_argument(
        "hmm",
        help="Path to HMM database file (for model sizes / work planning)",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output directory [default: {database}.BATHwater]",
    )
    parser.add_argument(
        "-p", "--processors",
        type=int,
        default=4,
        help="Controls partition granularity [default: 4]",
    )
    parser.add_argument(
        "--db-name",
        default=None,
        help="Filter to a specific database name in pass1_hits",
    )
    args = parser.parse_args()

    # Validate inputs
    if not os.path.isfile(args.database):
        log.error(f"Database not found: {args.database}")
        sys.exit(1)
    if not os.path.isfile(args.sequences):
        log.error(f"Sequences not found: {args.sequences}")
        sys.exit(1)
    if not os.path.isfile(args.hmm):
        log.error(f"HMM file not found: {args.hmm}")
        sys.exit(1)

    conn = sqlite3.connect(args.database)

    # Check pass1_hits exists and has data
    count = conn.execute("SELECT COUNT(*) FROM pass1_hits").fetchone()[0]
    if count == 0:
        log.error("No pass-1 hits in database. Run TEBinSorter first.")
        conn.close()
        sys.exit(1)
    log.info(f"Found {count} pass-1 hits in {args.database}")

    output_dir = args.output or f"{args.database}.BATHwater"

    emit_partitions(
        conn, args.sequences, args.hmm, output_dir,
        n_workers=args.processors,
        db_name=args.db_name,
    )

    conn.close()


if __name__ == "__main__":
    main()
