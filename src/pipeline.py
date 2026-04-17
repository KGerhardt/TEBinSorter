"""
Main pipeline for TE classification via two-pass HMM search.

Orchestrates: FASTA ingestion -> alphabet detection -> optional translation
-> pass 1 (coarse) -> optional pass 2 (sensitive) -> SQLite results.
"""

import argparse
import logging
import os
import time

from hmm import peek_alphabet, needs_translation, load_hmms, AMINO_ALPHABET, DNA_ALPHABET
from search import (build_sequence_block, pass1_screen, pass2_search,
                    pass2_search_parallel, legacy_search)
from sequence import translate_fasta, open_input
from results import (create_db, store_sequences, store_pass1, store_pass2,
                     store_legacy, export_tsv, export_best_hits_tsv,
                     export_all_domains_tsv, export_domain_sequences)
from emit import emit_partitions
from quick import quick_search
from iterative_search import iterative_search
from facet_classify import facet_classify, export_classifications_tsv
from id_registry import IDRegistry
from deconflict import (store_hits_numeric, load_hits_fast,
                        best_per_family_numeric, best_per_frame_numeric,
                        NUMERIC_HITS_SCHEMA)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Known database aliases -> paths relative to TEsorter database dir
DB_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..", "database",
)

DB_ALIASES = {
    "rexdb":    "REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm",
    "gydb":     "GyDB2.hmm",
    "line":     "Kapitonov_et_al.GENE.LINE.hmm",
    "tir":      "Yuan_and_Wessler.PNAS.TIR.hmm",
    "sine":     "AnnoSINE.hmm",
}


def resolve_db(name):
    """Resolve a database name or alias to an absolute path."""
    if os.path.isfile(name):
        return os.path.abspath(name)
    if name in DB_ALIASES:
        path = os.path.join(DB_DIR, DB_ALIASES[name])
        if os.path.isfile(path):
            return os.path.abspath(path)
        raise FileNotFoundError(f"Database alias '{name}' -> {path} not found")
    raise FileNotFoundError(f"Database '{name}' not found (not a file or known alias)")


def parse_args():
    parser = argparse.ArgumentParser(
        prog="ksort",
        description="TE classification via two-pass HMM search.",
    )
    parser.add_argument(
        "sequence",
        help="Input TE/LTR sequences in FASTA format",
    )
    parser.add_argument(
        "-d", "--database",
        default="rexdb",
        help="Comma-separated list of database names or paths "
             f"(aliases: {', '.join(DB_ALIASES.keys())}) [default: rexdb]",
    )
    parser.add_argument(
        "--max-search",
        action="store_true",
        default=False,
        help="Search against all known databases",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        default=False,
        help="Quick mode: sub-HMM triage assigns each frame to its best "
             "model, confirms with one nobias search, leftovers get full "
             "legacy search. Faster than default on typical inputs.",
    )
    parser.add_argument(
        "--facet",
        action="store_true",
        default=False,
        help="Facet mode: sub-HMM screen -> verify top hit per family -> "
             "confidence-tiered reporting -> legacy fallback. Fastest mode "
             "with ~94%% family recall.",
    )
    parser.add_argument(
        "--iterative",
        action="store_true",
        default=False,
        help="Iterative mode: facet screen -> graph-guided confirmation "
             "rounds -> legacy fallback on unresolved frames.",
    )
    parser.add_argument(
        "--two-pass",
        action="store_true",
        default=False,
        help="Use two-pass search: fast coarse screen with bias filter, "
             "then sensitive search on plausible pairs only. Faster on "
             "genomic input where most sequence is non-TE. Default mode "
             "is single-pass nobias (equivalent to TEsorter).",
    )
    parser.add_argument(
        "--pass-1-only",
        action="store_true",
        default=False,
        help="Run only the coarse pass-1 screen and store results. "
             "Implies --two-pass. Inspect the database before running pass 2.",
    )
    parser.add_argument(
        "-o", "--outdir",
        default=None,
        help="Output directory [default: {input}.TEBinSorter]",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output file prefix [default: basename of input]",
    )
    parser.add_argument(
        "-p", "--processors",
        type=int,
        default=4,
        help="Processors to use [default: 4]",
    )
    parser.add_argument(
        "--F1",
        type=float,
        default=0.02,
        help="MSV filter threshold for pass-1 screen. Higher values are "
             "more permissive. HMMER default is 0.02. Set to 0.1 to capture "
             "hits in compositionally biased sequences at ~10%% runtime cost. "
             "[default: %(default)s]",
    )
    parser.add_argument(
        "--emit-bath",
        action="store_true",
        default=False,
        help="After pass 1, emit routed FASTA partitions for the BATH "
             "aligner instead of running pass 2. Output goes to "
             "{prefix}.BATHwater/ directory.",
    )
    parser.add_argument(
        "--compat-tesorter-rounding",
        action="store_true",
        default=False,
        help="Round normalized scores to 2 decimal places before threshold "
             "comparison, replicating a TEsorter rounding bug. Use only for "
             "exact result reproduction against old TEsorter output.",
    )
    return parser.parse_args()


def run_database_legacy(db_path, seq_block, db_name, conn, registry=None):
    """
    Legacy mode: single-pass nobias search, all sequences against all models.
    """
    log.info(f"Loading HMMs from {db_name}")
    t0 = time.time()
    hmms = load_hmms(db_path)
    t1 = time.time()
    log.info(f"  Loaded {len(hmms)} models in {t1 - t0:.1f}s")

    if registry:
        registry.register_models(hmms)

    log.info(f"  Legacy search: bias filter OFF, all models, all sequences")
    t2 = time.time()
    hits = legacy_search(hmms, seq_block)
    t3 = time.time()
    log.info(f"  {len(hits)} hits in {t3 - t2:.1f}s")

    store_legacy(conn, hits, db_name, is_legacy=1)

    if registry:
        log.info(f"  Storing numeric hits")
        t4 = time.time()
        store_hits_numeric(conn, hits, registry)
        t5 = time.time()
        log.info(f"  Numeric storage in {t5 - t4:.1f}s")

    return len(hits)


def run_database(db_path, seq_block, seq_fasta, db_name, alphabet, conn,
                 pass1_only=False, n_workers=4, F1=0.02):
    """
    Run the two-pass search for a single database.

    Args:
        db_path: path to HMM database file
        seq_block: DigitalSequenceBlock (amino or nucl as appropriate)
        seq_fasta: path to the FASTA file (for parallel worker init)
        db_name: short name for this database (for tagging results)
        alphabet: easel.Alphabet for this database
        conn: sqlite3 connection for storing results
        pass1_only: if True, skip pass 2
        n_workers: number of worker processes for pass 2

    Returns:
        tuple of (pass1_hit_count, pass2_hit_count)
    """
    log.info(f"Loading HMMs from {db_name}")
    t0 = time.time()
    hmms = load_hmms(db_path)
    hmms_dict = {hmm.name: hmm for hmm in hmms}
    t1 = time.time()
    log.info(f"  Loaded {len(hmms)} models in {t1 - t0:.1f}s")

    # Pass 1
    log.info(f"  Pass 1: coarse screen (bias filter ON)")
    t2 = time.time()
    p1_hits, seq_models = pass1_screen(hmms, seq_block, F1=F1)
    t3 = time.time()
    n_seqs = len(seq_models)
    n_pairs = sum(len(v) for v in seq_models.values())
    log.info(f"  Pass 1: {len(p1_hits)} hits, {n_seqs} seqs with signal, "
             f"{n_pairs} seq-model pairs in {t3 - t2:.1f}s")

    store_pass1(conn, p1_hits, db_name)

    if pass1_only:
        log.info(f"  --pass-1-only: skipping pass 2 for {db_name}")
        return len(p1_hits), 0

    # Pass 2
    needed_models = set()
    for models in seq_models.values():
        needed_models |= models
    log.info(f"  Pass 2: sensitive search (bias filter OFF) on "
             f"{len(needed_models)} models")

    t4 = time.time()
    if n_workers > 1:
        p2_hits = pass2_search_parallel(
            db_path, seq_fasta, seq_models, hmms_dict,
            alphabet, n_workers=n_workers,
        )
    else:
        p2_hits = pass2_search(hmms_dict, seq_block, seq_models)
    t5 = time.time()
    log.info(f"  Pass 2: {len(p2_hits)} hits in {t5 - t4:.1f}s")

    store_pass2(conn, p2_hits, db_name)

    return len(p1_hits), len(p2_hits)


def main():
    args = parse_args()

    # Resolve output directory and prefix
    input_base = os.path.basename(args.sequence)
    outdir = args.outdir or f"{input_base}.TEBinSorter"
    prefix = args.prefix or input_base
    os.makedirs(outdir, exist_ok=True)

    db_path_out = os.path.join(outdir, f"{prefix}.db")
    aa_fasta = os.path.join(outdir, f"{prefix}.aa")

    # Resolve databases
    if args.max_search:
        db_names = list(DB_ALIASES.keys())
    else:
        db_names = [s.strip() for s in args.database.split(",")]

    db_paths = {}
    for name in db_names:
        db_paths[name] = resolve_db(name)

    log.info(f"Input: {args.sequence}")
    log.info(f"Databases: {', '.join(db_names)}")
    log.info(f"Output directory: {outdir}")
    log.info(f"File prefix: {prefix}")

    # Check which alphabets we need
    any_amino = False
    any_nucl = False
    db_alphabets = {}
    for name, path in db_paths.items():
        alphabet = peek_alphabet(path)
        db_alphabets[name] = alphabet
        if alphabet == AMINO_ALPHABET:
            any_amino = True
        else:
            any_nucl = True
        log.info(f"  {name}: {alphabet}, translate={'yes' if alphabet == AMINO_ALPHABET else 'no'}")

    # Create results database and ID registry
    conn = create_db(db_path_out)
    conn.executescript(NUMERIC_HITS_SCHEMA)
    registry = IDRegistry(conn)

    # Read input and store sequence metadata
    t_start = time.time()
    fa = open_input(args.sequence)
    nucl_lengths = {rec.name: len(rec) for rec in fa}
    store_sequences(conn, nucl_lengths)
    log.info(f"Input: {len(nucl_lengths)} sequences")

    # Translate if any database needs amino acid sequences
    aa_block = None
    if any_amino:
        # Remove stale index if present
        for f in [aa_fasta + ".fxi"]:
            if os.path.exists(f):
                os.remove(f)
        log.info("Six-frame translating input sequences")
        t0 = time.time()
        translate_fasta(args.sequence, aa_fasta)
        t1 = time.time()
        log.info(f"  Translation done in {t1 - t0:.1f}s -> {aa_fasta}")
        aa_block = build_sequence_block(aa_fasta, AMINO_ALPHABET)
        log.info(f"  Built amino acid sequence block: {len(aa_block)} frames")

    # Build nucleotide block if needed
    nucl_block = None
    if any_nucl:
        log.info("Building nucleotide sequence block")
        nucl_block = build_sequence_block(args.sequence, DNA_ALPHABET)
        log.info(f"  Built nucleotide sequence block: {len(nucl_block)} sequences")

    # Run each database
    for name in db_names:
        path = db_paths[name]
        alphabet = db_alphabets[name]

        if alphabet == AMINO_ALPHABET:
            seq_block = aa_block
            seq_fasta = aa_fasta
        else:
            seq_block = nucl_block
            seq_fasta = args.sequence

        log.info(f"--- Searching {name} ({os.path.basename(path)}) ---")

        two_pass = args.two_pass or args.pass_1_only or args.emit_bath

        if args.facet and alphabet != DNA_ALPHABET:
            t_f0 = time.time()
            classifications, f_legacy = facet_classify(
                path, seq_block, seq_fasta, alphabet,
                n_workers=args.processors,
                checkpoint_dir=outdir)
            t_f1 = time.time()
            n_primary = sum(1 for c in classifications if not c.get("is_secondary"))
            log.info(f"  Facet mode: {n_primary} assignments, "
                     f"{len(f_legacy)} legacy hits in {t_f1 - t_f0:.1f}s")
            # Store legacy fallback hits
            if f_legacy:
                store_legacy(conn, f_legacy, name, is_legacy=1)
            # Export classifications
            cls_tsv = os.path.join(outdir, f"{prefix}.{name}.classifications.tsv")
            export_classifications_tsv(classifications, cls_tsv)
            log.info(f"  Classifications: {cls_tsv}")
        elif args.facet and alphabet == DNA_ALPHABET:
            log.info(f"  DNA database: using legacy search (facets AA-only)")
            run_database_legacy(path, seq_block, name, conn, registry=registry)
        elif args.iterative:
            t_i0 = time.time()
            i_hits = iterative_search(
                path, seq_block, seq_fasta, alphabet,
                n_workers=args.processors,
                checkpoint_dir=outdir)
            t_i1 = time.time()
            log.info(f"  Iterative mode: {len(i_hits)} hits in {t_i1 - t_i0:.1f}s")
            store_legacy(conn, i_hits, name)
            if registry:
                store_hits_numeric(conn, i_hits, registry)
        elif args.quick:
            t_q0 = time.time()
            q_hits = quick_search(path, seq_block, seq_fasta, alphabet)
            t_q1 = time.time()
            log.info(f"  Quick mode: {len(q_hits)} hits in {t_q1 - t_q0:.1f}s")
            store_legacy(conn, q_hits, name)
        elif not two_pass:
            run_database_legacy(path, seq_block, name, conn, registry=registry)
        else:
            skip_pass2 = args.pass_1_only or args.emit_bath
            p1, p2 = run_database(
                path, seq_block, seq_fasta, name, alphabet, conn,
                pass1_only=skip_pass2,
                n_workers=args.processors,
                F1=args.F1,
            )

            if args.emit_bath:
                bath_dir = os.path.join(outdir, "BATHwater")
                log.info(f"  Emitting BATH partitions to {bath_dir}/")
                emit_partitions(conn, seq_fasta, path, bath_dir,
                                n_workers=args.processors, db_name=name)

    # TODO: rewrite exports with numpy for large datasets
    # Flat file exports temporarily disabled -- results are in the SQLite db
    # # Export flat files
    # log.info("Exporting results")
    #
    # # Determine which table to export from
    # if not two_pass:
    #     hit_table = "legacy_hits"
    # else:
    #     hit_table = "pass2_hits"
    #
    # def outpath(filename):
    #     return os.path.join(outdir, filename)
    #
    # if two_pass:
    #     p1_tsv = outpath(f"{prefix}.pass1.tsv")
    #     export_tsv(conn, p1_tsv, table="pass1_hits")
    #     log.info(f"  Raw pass-1 hits: {p1_tsv}")
    #
    # if not args.pass_1_only:
    #     raw_tsv = outpath(f"{prefix}.{'pass2' if two_pass else 'legacy'}.tsv")
    #     export_tsv(conn, raw_tsv, table=hit_table)
    #     log.info(f"  Raw hits: {raw_tsv}")
    #
    #     best_tsv = outpath(f"{prefix}.best.tsv")
    #     export_best_hits_tsv(conn, best_tsv, nucl_lengths=nucl_lengths,
    #                          table=hit_table)
    #     log.info(f"  Best hits: {best_tsv}")
    #
    #     domains_tsv = outpath(f"{prefix}.domains.tsv")
    #     export_all_domains_tsv(conn, domains_tsv, nucl_lengths=nucl_lengths,
    #                            table=hit_table)
    #     log.info(f"  All domains: {domains_tsv}")
    #
    #     if any_amino:
    #         dom_faa = outpath(f"{prefix}.domains.faa")
    #         export_domain_sequences(conn, dom_faa, aa_fasta,
    #                                 nucl_lengths=nucl_lengths,
    #                                 table=hit_table)
    #         log.info(f"  Domain sequences: {dom_faa}")

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results database: {db_path_out}")

    conn.close()


if __name__ == "__main__":
    main()
