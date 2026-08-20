"""
Main pipeline for TE classification.

Orchestrates: FASTA ingestion -> alphabet detection -> optional translation
-> HMM search -> classification -> BLAST pass-2 -> SQLite + TSV output.
"""

import argparse
import logging
import os
import sys
import time

from .paths import get_db_dir
from .hmm import peek_alphabet, needs_translation, load_hmms, AMINO_ALPHABET, DNA_ALPHABET
from .search import build_sequence_block, legacy_search, legacy_search_nucl
from .sequence import translate_fasta, open_input
from .results import (create_db, store_sequences, store_legacy,
                     index_hits_tables, finalize_db)
from .classifier import (classify_sequences, export_classification_tsv,
                       store_classifications, reconcile_classifications,
                       DB_CONFIGS)
from .blast_pass2 import blast_pass2
from .hierarchical_search import ENGINES
from . import bath_search


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
log = logging.getLogger(__name__)

# Database directory: the databases ship inside the package. Override with
# --db-dir or TESORTER2_DB to point at a custom collection.
DB_DIR = get_db_dir()

DB_ALIASES = {
    "rexdb":    "REXdb_protein_database_viridiplantae_v4.0_plus_metazoa_v3.1.hmm",
    "gydb":     "GyDB2.hmm",
    "line":     "Kapitonov_et_al.GENE.LINE.hmm",
    "tir":      "Yuan_and_Wessler.PNAS.TIR.hmm",
    "sine":     "AnnoSINE_core.hmm",
    "sine-so":  "SINE_SO.hmm",
}


def resolve_db(name, db_dir=None):
    """Resolve a database name or alias to an absolute path."""
    if os.path.isfile(name):
        return os.path.abspath(name)
    if name in DB_ALIASES:
        base = db_dir or DB_DIR
        path = os.path.join(base, DB_ALIASES[name])
        if os.path.isfile(path):
            return os.path.abspath(path)
        raise FileNotFoundError(
            f"Database alias '{name}' -> {path} not found under {base}")
    raise FileNotFoundError(f"Database '{name}' not found (not a file or known alias)")


# Subcommands. The mode is the verb: what the input *is* decides almost
# everything downstream -- which engines can run, whether the input is
# windowed, whether there is a per-element classification at all -- so it reads
# better as a mode than as a flag hung off a single command.
MODES = ("sequences", "genome")


def _add_shared(p):
    """Options meaningful in both modes."""
    p.add_argument("sequence", help="Input FASTA")
    p.add_argument(
        "-d", "--database", default=None,
        help="Comma-separated list of database names or paths "
             f"(aliases: {', '.join(DB_ALIASES.keys())}). "
             "[default: every bundled database except sine-so]")
    p.add_argument("--max-search", action="store_true", default=False,
                   help="Search against all known databases. This is the "
                        "default when -d is omitted; the flag still forces the "
                        "full set alongside an explicit -d.")
    p.add_argument("-o", "--outdir", default=None,
                   help="Output directory [default: {input}.TEsorter2]")
    p.add_argument("--prefix", default=None,
                   help="Output file prefix [default: basename of input]")
    p.add_argument("--db-dir", default=None,
                   help="Directory containing the HMM databases "
                        "[default: the databases bundled with the package]")
    p.add_argument("-p", "--processors", type=int, default=4,
                   help="Processors to use [default: 4]")
    p.add_argument("--include-sine-so", action="store_true", default=False,
                   help="Include the SINE_SO model (M=4176) in AnnoSINE "
                        "searches. Excluded by default: it costs 71%% of "
                        "AnnoSINE's compute for <1 filterable hit per 100k "
                        "sequences.")
    p.add_argument("--mask-stops", action="store_true", default=False,
                   help="Six-frame translate stop codons as 'X' (unknown "
                        "residue) rather than '*' (terminator), letting a "
                        "profile align through a premature stop instead of "
                        "being cut short by it. Recovers domains in degraded "
                        "copies, but covers only the in-frame part of what "
                        "BATH does and its false-positive cost is not yet "
                        "characterized.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="tesorter2",
        description="Fast TE classification using HMM profile databases.")
    sub = parser.add_subparsers(dest="mode", metavar="{sequences,genome}")

    # ---- tesorter2 sequences ----
    seq = sub.add_parser(
        "sequences", help="Classify pre-extracted TE sequences",
        description="Classify pre-extracted TE/LTR sequences. Emits a "
                    "per-element classification per database, a reconciled "
                    "combined call, and a BLAST pass-2 over whatever no "
                    "database classified.")
    _add_shared(seq)
    seq.add_argument(
        "--stages", default=None,
        help="Run the staged search cascade over amino-acid databases instead "
             "of a single engine: an ordered, comma-separated engine list "
             f"({', '.join(sorted(k for k in ENGINES if k != 'nhmmer'))}). "
             "Each stage searches only what earlier stages left unresolved, "
             "per database. DNA databases always use nhmmer, the only engine "
             "that reads DNA profiles. e.g. --stages nail,hmmer,bath. "
             "A cascade containing nail also requires --mask-stops.")
    seq.add_argument("--bath", action="store_true", default=False,
                     help="Use BATH (frameshift-aware translated nucleotide "
                          "search) instead of hmmer for amino-acid databases. "
                          "Set BATH_BIN_DIR if bathsearch/bathconvert are not "
                          "on PATH.")
    seq.add_argument("--dna-engine", choices=("nhmmer", "hmmsearch"),
                     default="nhmmer",
                     help="Engine for DNA databases. nhmmer scans both strands "
                          "and handles long targets; hmmsearch scores only the "
                          "strand it is given. [default: nhmmer]")
    seq.add_argument("--compat-tesorter-voting", action="store_true",
                     default=False,
                     help="Decide the clade by raw domain-count plurality, "
                          "replicating TEsorter exactly. Default is a "
                          "score-weighted clade vote.")
    seq.add_argument("--compat-tesorter-rounding", action="store_true",
                     default=False,
                     help="Round normalized scores to 2 dp before threshold "
                          "comparison, replicating a TEsorter rounding bug.")
    seq.add_argument("--compat-tesorter-output", action="store_true",
                     default=False,
                     help="Emit the combined .cls.tsv in TEsorter's 7-column "
                          "format, without SecondaryHits/SO/lineage.")
    seq.add_argument("--no-tesorter-outputs", action="store_true",
                     default=False,
                     help="Skip the TEsorter-compatible companion files "
                          "(.cls.lib, .cls.pep, .dom.gff3, .dom.tsv, .dom.faa).")
    seq.add_argument("-nolib", "--no-library", action="store_true",
                     default=False,
                     help="Do not write the RepeatMasker library (.cls.lib).")
    seq.add_argument("-norc", "--no-reverse", action="store_true",
                     default=False,
                     help="Do not reverse-complement minus-strand sequences "
                          "when writing the .cls.lib library.")
    seq.add_argument(
        "--min-clade-delta", type=float, default=0.0, metavar="BITS",
        help="Demote a clade call to 'mixture' when the vote could be "
             "overturned by a per-hit score offset smaller than BITS. The "
             "clade vote sums score/model_len, so a uniform offset favours "
             "clades carried by short models; this is the margin the call can "
             "absorb. Only the Clade field is affected -- Order and "
             "Superfamily are never touched. 0 (the default) never demotes, "
             "leaving classification "
             "bit-for-bit unchanged. ~5.6 is the measured pyhmmer-to-BATH "
             "offset, so that is the scale at which a call is engine-dependent.")
    seq.add_argument(
        "--emit-clade-delta", action="store_true", default=False,
        help="Append a CladeDelta column to the per-database .cls.tsv with "
             "that margin in bits ('inf' = only one clade in play). "
             "Diagnostic only; does not change any call.")
    seq.add_argument("--emit-bath", action="store_true", default=False,
                     help="Emit routed FASTA partitions for BATH to "
                          "{outdir}/BATHwater/.")
    # Deprecated, hidden
    for dep in ("--quick", "--iterative", "--two-pass", "--pass-1-only"):
        seq.add_argument(dep, action="store_true", default=False,
                         help=argparse.SUPPRESS)
    seq.add_argument("--F1", type=float, default=0.02, help=argparse.SUPPRESS)

    # ---- tesorter2 genome ----
    gen = sub.add_parser(
        "genome", help="Annotate TE domains across a whole assembly",
        description="Scan whole-genome sequence for TE domains. Windows the "
                    "genome, classifies each domain occurrence individually, "
                    "and emits a feature-level GFF3 + summary. No per-element "
                    "classification and no BLAST pass-2. Protein profiles are "
                    "searched with BATH (required); DNA profiles with stock "
                    "HMMER nhmmer.")
    _add_shared(gen)
    gen.add_argument("--win-size", type=float, default=1e6,
                     help="Window size for chunking [default: 1e6]")
    gen.add_argument("--win-ovl", type=float, default=1e5,
                     help="Window overlap [default: 1e5]")

    args = parser.parse_args(_normalize_argv(argv))
    if args.mode is None:
        parser.error("a mode is required: %s" % " | ".join(MODES))
    return args


def _normalize_argv(argv):
    """Accept the pre-subcommand form, so existing scripts keep working.

    `tesorter2 in.fa --genome ...` becomes `tesorter2 genome in.fa ...`, and a
    bare `tesorter2 in.fa ...` becomes `tesorter2 sequences in.fa ...`.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in MODES or argv[0] in ("-h", "--help"):
        return argv
    if argv[0].startswith("-"):
        return argv                      # let argparse report it
    if "--genome" in argv:
        log.warning("`--genome` is deprecated; use `tesorter2 genome ...`")
        return ["genome"] + [a for a in argv if a != "--genome"]
    log.warning("Assuming `tesorter2 sequences ...`; name the mode explicitly.")
    return ["sequences"] + argv


def run_database_legacy(db_path, seq_block, db_name, conn, alphabet=None,
                        dna_engine="nhmmer", n_workers=0):
    """
    Exhaustive single-pass nobias search against all models.

    DNA-alphabet databases go through nhmmer, which scans both strands;
    amino-acid databases go through hmmsearch. Pass dna_engine="hmmsearch" to
    force the old single-strand behaviour on DNA databases.

    Hits go to legacy_hits, which is now the only hits table.
    """
    log.info(f"Loading HMMs from {db_name}")
    t0 = time.time()
    hmms = load_hmms(db_path)
    use_nhmmer = alphabet == DNA_ALPHABET and dna_engine == "nhmmer"
    optimized = None
    if not use_nhmmer:
        # nhmmer's pipeline builds its own profiles; optimizing here would be
        # wasted work.
        from .hmm import build_optimized_profiles
        optimized = build_optimized_profiles(hmms, alphabet=alphabet)
    t1 = time.time()
    log.info(f"  Loaded {len(hmms)} models in {t1 - t0:.1f}s")

    t2 = time.time()
    if use_nhmmer:
        log.info(f"  nhmmer search: bias filter OFF, both strands, all models")
        hits = legacy_search_nucl(hmms, seq_block, n_workers=n_workers)
    else:
        log.info(f"  Legacy search: bias filter OFF, all models, all sequences")
        hits = legacy_search(hmms, seq_block, optimized=optimized)
    t3 = time.time()
    log.info(f"  {len(hits)} hits in {t3 - t2:.1f}s")

    store_legacy(conn, hits, db_name)
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

    genome_mode = args.mode == "genome"

    if not genome_mode and args.stages:
        if args.two_pass or args.pass_1_only:
            raise SystemExit(
                "--stages is incompatible with the deprecated two-pass modes.")

    # Resolve output directory and prefix
    input_base = os.path.basename(args.sequence)
    outdir = args.outdir or f"{input_base}.TEsorter2"
    prefix = args.prefix or input_base
    os.makedirs(outdir, exist_ok=True)

    db_path_out = os.path.join(outdir, f"{prefix}.db")
    aa_fasta = os.path.join(outdir, f"{prefix}.aa")

    # Resolve databases. Omitting -d searches everything: each database
    # classifies independently and reconciliation resolves them afterwards, so
    # more databases is more evidence rather than more ambiguity.
    if args.max_search or not args.database:
        db_names = [k for k in DB_ALIASES.keys() if k != "sine-so"]
        if args.include_sine_so:
            db_names.append("sine-so")
    else:
        db_names = [s.strip() for s in args.database.split(",")]

    db_dir = get_db_dir(args.db_dir)
    db_paths = {}
    for name in db_names:
        db_paths[name] = resolve_db(name, db_dir=db_dir)

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

    # Genome mode: dispatch to the domain-level genome scanner and stop here.
    # It windows the genome, finds/classifies each TE protein domain, and emits
    # a GFF3 + summary -- no element classification, reconcile, or BLAST pass-2.
    if genome_mode:
        from . import genome

        log.info("Genome mode")
        genome.run_genome(
            args.sequence, db_paths, db_alphabets, outdir, prefix,
            n_workers=args.processors,
            win_size=args.win_size, win_ovl=args.win_ovl)
        return

    # Create results database
    conn = create_db(db_path_out)

    # Read input and store sequence metadata
    t_start = time.time()
    fa = open_input(args.sequence)
    nucl_lengths = {rec.name: len(rec) for rec in fa}
    store_sequences(conn, nucl_lengths)
    log.info(f"Input: {len(nucl_lengths)} sequences")

    # Translate if any database needs amino acid sequences. Under --bath the
    # amino-acid databases are searched by bathsearch directly against the
    # nucleotide input (frameshift-aware translated search), so no six-frame
    # translation is needed for them.
    aa_block = None
    if any_amino and not args.bath and not args.stages:
        # Remove stale index if present
        for f in [aa_fasta + ".fxi"]:
            if os.path.exists(f):
                os.remove(f)
        log.info("Six-frame translating input sequences%s",
                 " (stops masked as X)" if args.mask_stops else "")
        t0 = time.time()
        translate_fasta(args.sequence, aa_fasta, mask_stops=args.mask_stops)
        t1 = time.time()
        log.info(f"  Translation done in {t1 - t0:.1f}s -> {aa_fasta}")
        aa_block = build_sequence_block(aa_fasta, AMINO_ALPHABET)
        log.info(f"  Built amino acid sequence block: {len(aa_block)} frames")

    # Build nucleotide block if needed
    nucl_block = None
    if any_nucl and not args.stages:
        log.info("Building nucleotide sequence block")
        nucl_block = build_sequence_block(args.sequence, DNA_ALPHABET)
        log.info(f"  Built nucleotide sequence block: {len(nucl_block)} sequences")

    per_db_results = {}

    # --- Staged cascade (--stages) ---
    # Replaces both the per-database search loop and the classification loop
    # below; everything after them -- cross-database reconciliation, BLAST
    # pass-2, the combined export -- is shared and runs unchanged.
    if args.stages:
        from . import hierarchical_search

        per_db_results = hierarchical_search.run_cascade(
            conn, args.sequence, db_paths, db_alphabets, outdir,
            protein_stages=[x.strip() for x in args.stages.split(",")],
            n_workers=args.processors,
            compat_rounding=args.compat_tesorter_rounding,
            compat_voting=args.compat_tesorter_voting,
            mask_stops=args.mask_stops,
            min_clade_delta=getattr(args, "min_clade_delta", 0.0))

        log.info("Indexing hits tables")
        index_hits_tables(conn)

        for name, results in per_db_results.items():
            if not results:
                continue
            cls_tsv = os.path.join(outdir, f"{prefix}.{name}.cls.tsv")
            export_classification_tsv(
                results, cls_tsv,
                include_clade_delta=getattr(args, "emit_clade_delta", False))
            log.info(f"  {name}: {len(results)} classified -> {cls_tsv}")

            # Companion files: the RepeatMasker library only. The domain-level
            # files (.dom.gff3/.dom.tsv/.dom.faa, .cls.pep) encode coordinates
            # on the six-frame translation, and a cascade's hits for one
            # database can come from several engines -- BATH reports nucleotide
            # coordinates, nail its own -- so there is no single frame to write
            # them against. Same restriction --bath already carries.
            if not args.no_tesorter_outputs:
                from .tesorter_output import generate_all_outputs
                generate_all_outputs(
                    conn, os.path.join(outdir, f"{prefix}.{name}"), name,
                    args.sequence, None, nucl_lengths, results,
                    seq_type="nucl",
                    no_reverse=args.no_reverse,
                    no_library=args.no_library,
                    hits_table="legacy_hits",
                    domain_files=False,
                    keep=None,
                )

    # Databases handled by the single-engine path. Empty under --stages, which
    # has already classified everything through the cascade above.
    single_engine_dbs = [] if args.stages else db_names

    # Run each database
    for name in single_engine_dbs:
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

        if args.bath and alphabet == AMINO_ALPHABET:
            t_b0 = time.time()
            hits = bath_search.run_and_parse(
                path, args.sequence, name,
                n_workers=args.processors, outdir=outdir)
            store_legacy(conn, hits, name)
            log.info(f"  BATH search: {len(hits)} hits in {time.time() - t_b0:.1f}s")
        else:
            run_database_legacy(path, seq_block, name, conn, alphabet=alphabet,
                                dna_engine=args.dna_engine,
                                n_workers=args.processors)

    # Build hits-table indexes now that all HMM hits are written and
    # before the classification phase starts reading from them.
    log.info("Indexing hits tables")
    t_idx0 = time.time()
    index_hits_tables(conn)
    log.info(f"  Indexed in {time.time() - t_idx0:.1f}s")

    # --- Classification ---
    log.info("--- Classification ---")
    from .deconflict import load_hits

    for name in single_engine_dbs:
        config = DB_CONFIGS.get(name)
        if config is None:
            log.warning(f"  No classifier config for {name}, skipping")
            continue

        hits_table = "legacy_hits"
        hits = load_hits(db_path_out, table=hits_table, database=name)
        if hits is None:
            continue

        log.info(f"  Classifying {name} ({hits_table})")
        results = classify_sequences(hits, config,
                                     compat_rounding=args.compat_tesorter_rounding,
                                     compat_voting=args.compat_tesorter_voting,
                                     min_clade_delta=getattr(args, "min_clade_delta", 0.0))

        # Store and export per-database classification (TEsorter format)
        store_classifications(conn, results, database=name)
        cls_tsv = os.path.join(outdir, f"{prefix}.{name}.cls.tsv")
        export_classification_tsv(
            results, cls_tsv,
            include_clade_delta=getattr(args, "emit_clade_delta", False))
        log.info(f"    {len(results)} classified -> {cls_tsv}")

        # TEsorter-compatible companion files, named {prefix}.{db}.* to sit
        # beside the per-database .cls.tsv. The domain-level files encode
        # six-frame translated coordinates, so they are only written when this
        # database was searched on the translated block; --bath and
        # DNA-alphabet databases get the library alone.
        if not args.no_tesorter_outputs and results:
            from .tesorter_output import generate_all_outputs, domain_keys
            from .classifier import select_domain_indices
            six_frame = (aa_block is not None
                         and db_alphabets[name] == AMINO_ALPHABET)
            keep = domain_keys(hits, select_domain_indices(
                hits, config, compat_rounding=args.compat_tesorter_rounding))
            generate_all_outputs(
                conn, os.path.join(outdir, f"{prefix}.{name}"), name,
                args.sequence,
                aa_fasta if six_frame else None,
                nucl_lengths, results,
                seq_type="nucl",
                no_reverse=args.no_reverse,
                no_library=args.no_library,
                hits_table=hits_table,
                domain_files=six_frame,
                keep=keep,
            )

        per_db_results[name] = results

    # Reconcile across databases via hierarchical weighted vote
    reconciled = reconcile_classifications(per_db_results)
    all_classifications = {r["id"]: r for r in reconciled}
    log.info(f"  Reconciled across {len(per_db_results)} databases: "
             f"{len(reconciled)} sequences")

    # --- BLAST pass-2 ---
    all_results = list(reconciled)
    if not args.pass_1_only and all_classifications:
        log.info("--- BLAST pass-2 ---")
        blast_cls = blast_pass2(
            args.sequence, conn,
            hmm_classifications=all_classifications,
            seq_type="nucl",
            n_processors=args.processors,
            outdir=outdir,
        )

        if blast_cls:
            store_classifications(conn, blast_cls, database="blast_pass2")
            all_results = list(reconciled) + blast_cls

    # Export combined classification
    if all_results:
        combined_tsv = os.path.join(outdir, f"{prefix}.cls.tsv")
        # Lineage rides on the combined file only: the per-database .cls.tsv
        # stays at TEsorter's exact 7 columns, which is what EDTA consumes.
        export_classification_tsv(
            all_results, combined_tsv,
            include_secondary=not args.compat_tesorter_output,
            include_so=not args.compat_tesorter_output,
            include_lineage=not args.compat_tesorter_output,
        )
        log.info(f"  Combined: {len(all_results)} classified -> {combined_tsv}")

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

    # Build remaining indexes (classifications, blast_hits) and truncate WAL
    log.info("Finalizing database")
    t_fin0 = time.time()
    finalize_db(conn)
    log.info(f"  Finalized in {time.time() - t_fin0:.1f}s")

    t_end = time.time()
    log.info(f"Done in {t_end - t_start:.1f}s")
    log.info(f"Results database: {db_path_out}")

    conn.close()


if __name__ == "__main__":
    main()
