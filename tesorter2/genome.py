"""
genome.py — Genome mode: scan whole-genome sequences for TE protein domains.

Mirrors TEsorter's `-genome` mode (TEsorter/app.py: genomeAnn / cut_seqs /
_hmm2best(genome=True) / hmm2best(genome=True) / group_resolve_overlaps /
summary_genome). The genome is windowed into overlapping fragments; each window
is searched for TE protein-domain models; EACH domain occurrence is classified
individually; coordinates are mapped back to genomic space; overlapping features
are resolved; and a domain-level GFF3 + summary table are emitted.

Unlike element mode (classifier.classify_sequences), genome mode:
  * does NOT collapse domains into one classification per input sequence —
    every domain occurrence becomes its own GFF feature,
  * runs NO BLAST pass-2,
  * emits NO per-element .cls.tsv.

Two engines, selected by --bath:
  * HMMER (default): six-frame translate each window (sequence.translate_fasta),
    legacy_search, then map amino-acid envelope coords -> nucleotide coords via
    sequence.parse_frame_suffix + sequence.aa_to_nucl_coords.
  * BATH (--bath): bathsearch --fs directly on the nucleotide windows; the
    tblout already yields nucleotide ali coords + strand per domain occurrence
    (ascending-normalized, |fwd1/|rev1 suffix), so no translation / frame math.
"""

import itertools
import logging
import os
import re
import sys
import time

from .hmm import load_hmms, build_optimized_profiles, AMINO_ALPHABET
from .sequence import (open_input, clean_seq, revcomp, translate_fasta,
                      parse_frame_suffix, aa_to_nucl_coords, load_sequences_dict)
from .search import build_sequence_block, legacy_search
from .so_map import so_gff_type
from .classifier import (parse_clade_rexdb, parse_clade_gydb, classify_element,
                        load_gydb_clade_map, DB_CONFIGS)
from . import bath_search


log = logging.getLogger(__name__)

# Order of TE orders for summary sorting (mirrors TEsorter app.py:56).
ORDERS = ['LTR', 'pararetrovirus', 'DIRS', 'Penelope', 'LINE', 'SINE',
          'TIR', 'Helitron', 'Maverick', 'mixture', 'Unknown', 'Total']

# Hit filters — match classifier.apply_filters defaults / TEsorter genome mode.
MIN_COV = 20.0
MAX_EVALUE = 1e-3
MIN_ACC = 0.5
MIN_NORM_SCORE = 0.1

# pyhmmer's search pipeline rejects single sequences longer than 100k residues.
# TEsorter avoids this by calling the hmmscan binary; our HMMER engine uses
# pyhmmer, so we cap each window so its translated frame stays under the limit.
# The cap is invisible to results: the retained overlap (>= any TE ORF) means
# every domain is still fully contained in at least one window. BATH has no such
# limit (it searches nucleotides directly) and keeps the user's window size.
_PYHMMER_MAX_AA = 100000
_HMMER_MAX_WIN_NUCL = 240000   # ~80k AA per frame, comfortably under the limit
_HMMER_MAX_WIN_OVL = 60000     # >> longest TE protein domain

# Window IDs are "{chrom}:{start}-{end}" with a 1-based start (see cut_genome).
# rsplit on the LAST ':' so chromosome names containing ':' survive.
_WIN_RE = re.compile(r'^(\d+)-(\d+)$')


# ---------------------------------------------------------------------------
# Windowing (ports TEsorter modules/split_records.cut_seqs)
# ---------------------------------------------------------------------------

def cut_genome(in_fasta, out_fasta, win_size=int(1e6), win_ovl=int(1e5)):
    """Window genome sequences into overlapping fragments.

    Step = win_size; each window spans win_size + win_ovl. Window FASTA IDs are
    '{chrom}:{start}-{end}' with a 1-based start, matching TEsorter so the
    coordinate offset can be parsed straight back out.

    Returns {window_id: window_nucleotide_length}.
    """
    win_size, win_ovl = int(win_size), int(win_ovl)
    win_lengths = {}
    # Drop any stale pyfastx index so downstream readers don't see old windows.
    if os.path.exists(out_fasta + ".fxi"):
        os.remove(out_fasta + ".fxi")
    with open(out_fasta, "w") as fout:
        for rec in open_input(in_fasta):
            seq = clean_seq(str(rec.seq).upper())
            seq_len = len(seq)
            for s in range(0, seq_len + 1, win_size):
                e = min(s + win_size + win_ovl, seq_len)
                sseq = seq[s:e]
                if not sseq:
                    continue
                wid = "{}:{}-{}".format(rec.name, s + 1, e)
                win_lengths[wid] = len(sseq)
                fout.write(">{}\n{}\n".format(wid, sseq))
                if e == seq_len:
                    break
    return win_lengths


def _split_window_id(window_id):
    """('chrom:1-1000000') -> ('chrom', 1, 1000000). Offset is the 1-based
    genomic start of the window. Returns (window_id, None, None) if not a
    windowed id (so non-windowed input degrades gracefully)."""
    base, _, span = window_id.rpartition(":")
    if base:
        m = _WIN_RE.match(span)
        if m:
            return base, int(m.group(1)), int(m.group(2))
    return window_id, None, None


# ---------------------------------------------------------------------------
# Overlap resolution (ports TEsorter app.py:843-885)
# ---------------------------------------------------------------------------

def _overlap_pct(a, b):
    """Percent overlap relative to the shorter feature. Features are tuples
    with start at index 3 and end at index 4."""
    ovl = max(0, min(a[4], b[4]) - max(a[3], b[3]))
    shorter = min((a[4] - a[3] + 1), (b[4] - b[3] + 1))
    return 100.0 * ovl / shorter if shorter > 0 else 0.0


def resolve_overlaps(features, max_ovl=20):
    """Within one chromosome, drop equal/overlapping features keeping the
    higher score (feature[5]). Faithful port of TEsorter resolve_overlaps."""
    last = None
    discards = []
    for feat in sorted(features, key=lambda x: x[3]):
        discard = None
        if last:
            if feat == last:
                pair = [last, feat]            # retain, discard (duplicate)
            elif _overlap_pct(feat, last) > max_ovl:
                pair = [feat, last] if feat[5] > last[5] else [last, feat]
            else:
                last = feat
                continue
            _, discard = pair
            discards.append(discard)
        if not last or discard != feat:
            last = feat
    return sorted(set(features) - set(discards), key=lambda x: x[3])


def group_resolve_overlaps(features):
    """Resolve overlaps per chromosome (feature[0])."""
    resolved = []
    for chrom, items in itertools.groupby(sorted(features, key=lambda x: x[0]),
                                           key=lambda x: x[0]):
        items = list(items)
        log.info("  resolving overlaps in {} ({} features)".format(chrom, len(items)))
        resolved += resolve_overlaps(items)
    return resolved


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def _fmt_cls(*args):
    """Join classification levels, dropping 'unknown' and duplicates
    (faithful port of TEsorter fmt_cls)."""
    values = []
    for arg in args:
        if arg == "unknown" or arg in values:
            continue
        values.append(arg)
    return "/".join(values)


def _parse_model(model, config):
    """Model name -> (gene, clade) for the configured parser."""
    parser = config["clade_parser"]
    if parser == "rexdb":
        _, _, clade, gene = parse_clade_rexdb(model)
        return gene, clade
    if parser == "gydb":
        return parse_clade_gydb(model)
    return "SINE", "SINE"


# ---------------------------------------------------------------------------
# Per-hit -> genomic feature
# ---------------------------------------------------------------------------

def _hit_to_feature(h, engine, config, win_lengths):
    """Turn one passing domain hit into a GFF feature tuple, classified on its
    own. Returns (feature_tuple, aa_or_nucl_domain_seq_key) or None.

    Feature tuple layout (indices used by resolve_overlaps): 0 chrom, 1 source,
    2 type, 3 start, 4 end, 5 score, 6 strand, 7 frame, 8 attributes,
    9 gid, 10 (target_name, w_start, w_end) for sequence extraction.
    """
    model = h["query_name"]
    model_len = h["query_len"]
    if model_len <= 0:
        return None
    # Per-DOMAIN score (TEsorter genome uses domscore/tlen). A genome window
    # holds many domains per model, so the sequence-level score (col 7) would
    # be inflated; dom_score is this domain alone. For BATH, dom_score==score.
    norm_score = h["dom_score"] / model_len
    hmm_cov = 100.0 * (h["hmm_to"] - h["hmm_from"] + 1) / model_len
    if not (hmm_cov >= MIN_COV and h["evalue"] <= MAX_EVALUE
            and h["acc"] >= MIN_ACC and norm_score >= MIN_NORM_SCORE):
        return None

    gene, clade = _parse_model(model, config)

    if engine == "bath":
        window_id, strand, _frame = parse_frame_suffix(h["target_name"])
        w_start, w_end = h["ali_from"], h["ali_to"]   # nucleotide, ascending
        frame_str = "."
        # Domain sequence is the nucleotide window slice (revcomp on minus).
        extract = (window_id, w_start, w_end, strand)
    else:  # hmmer: AA envelope -> within-window nucleotide coords
        window_id, strand, frame = parse_frame_suffix(h["target_name"])
        wlen = win_lengths.get(window_id)
        if wlen is None:
            return None
        w_start, w_end = aa_to_nucl_coords(
            h["env_from"], h["env_to"], strand, frame, wlen)
        frame_str = str(frame)
        # Domain sequence is the AA slice of the translated frame.
        extract = (h["target_name"], h["env_from"], h["env_to"], "+")

    chrom, offset, _ = _split_window_id(window_id)
    if offset is not None:
        g_start = offset + w_start - 1
        g_end = offset + w_end - 1
    else:
        g_start, g_end = w_start, w_end
    if g_start > g_end:
        g_start, g_end = g_end, g_start

    order, superfamily, disp_clade, complete = classify_element(
        [gene], [clade], [model], [norm_score], config)
    cls = _fmt_cls(order, superfamily, disp_clade)

    gid = "{}:{}-{}|{}".format(chrom, g_start, g_end, model)
    name = "{}-{}".format(clade, gene)
    # The feature is a protein domain, not the element, so its type stays CDS:
    # typing it as the element's SO term would assert the domain *is* the
    # retrotransposon. The element's ontology term goes in Ontology_term, the
    # GFF3-reserved attribute for ontology cross-references.
    so_name, so_id = so_gff_type(order, superfamily)
    attr = ("ID={};Name={};Classification={};gene={};clade={};"
            "coverage={:.1f};evalue={:g};score={:.3f};"
            "Ontology_term={};so_name={}").format(
        gid, name, cls, gene, clade, hmm_cov, h["evalue"], norm_score,
        so_id, so_name)

    feature = (chrom, "TEsorter2", "CDS", g_start, g_end,
               round(norm_score, 3), strand, frame_str, attr,
               gid, extract)
    return feature, cls


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _extract_domain_seq(feature, engine, seq_cache):
    """Domain sub-sequence: amino acids for HMMER (AA slice of the translated
    frame), nucleotides for BATH (window slice, revcomp on minus). The
    extraction tuple is (key, start, end, strand) with key into seq_cache and
    1-based start/end in that sequence's own coordinate system."""
    key, start, end, strand = feature[10]
    seq = seq_cache.get(key, "")
    sub = seq[start - 1:end]
    if strand == "-":
        sub = revcomp(sub)
    return sub


def write_outputs(features, engine, outdir, prefix, seq_path):
    """Write {prefix}.dom.gff3, {prefix}.dom.{faa|fna}, and the summary table.
    Returns (gff_path, seq_out_path, summary_path)."""
    gff_path = os.path.join(outdir, "{}.dom.gff3".format(prefix))
    seq_ext = "faa" if engine == "hmmer" else "fna"
    seq_out_path = os.path.join(outdir, "{}.dom.{}".format(prefix, seq_ext))
    summary_path = os.path.join(outdir, "{}.genome.summary.tsv".format(prefix))

    # For sequence extraction: HMMER needs the AA frames, BATH the nucl windows.
    seq_cache = load_sequences_dict(seq_path)

    with open(gff_path, "w") as fg, open(seq_out_path, "w") as fs:
        fg.write("##gff-version 3\n")
        for feat in features:
            fg.write("\t".join(map(str, feat[:9])) + "\n")
            sub = _extract_domain_seq(feat, engine, seq_cache)
            fs.write(">{} {}\n{}\n".format(feat[9], feat[8], sub))

    summarize(features, summary_path)
    return gff_path, seq_out_path, summary_path


def summarize(features, summary_path):
    """Collapse classified features into a per-classification stats table
    (Order/Superfamily/Clade/Number/Total_length/Mean_length). Ports the
    spirit of TEsorter summary_genome, robust to orders outside ORDERS."""
    d_stats = {}
    for feat in sorted(features, key=lambda x: (x[0], x[3])):
        # Classification is stored in the attribute string.
        m = re.search(r'Classification=([^;]+)', feat[8])
        if not m:
            continue
        cls = tuple(m.group(1).split("/"))
        length = feat[4] - feat[3] + 1
        # expand to Order, Order/SF, Order/SF/Clade, and Total
        xcls = []
        if len(cls) == 3:
            xcls += [cls[:1], cls[:2]]
        elif len(cls) == 2:
            xcls += [cls[:1]]
        xcls += [cls, ("Total",)]
        for c in xcls:
            d_stats.setdefault(c, []).append(length)

    def order_key(item):
        c = item[0]
        try:
            oi = ORDERS.index(c[0])
        except ValueError:
            oi = len(ORDERS)
        return (oi, c)

    lines = [["Order", "Superfamily", "Clade", "Number",
              "Total_length", "Mean_length"]]
    for cls, lens in sorted(d_stats.items(), key=order_key):
        padded = [""] * (len(cls) - 1) + [cls[-1]] + [""] * (3 - len(cls))
        n, total = len(lens), sum(lens)
        lines.append(padded + [str(n), str(total), str(round(total / n, 1))])

    with open(summary_path, "w") as f:
        for line in lines:
            f.write("\t".join(line) + "\n")

    # Echo to stdout (TEsorter prints the summary).
    log.info("Summary of classifications:")
    for line in lines:
        sys.stdout.write("\t".join(line) + "\n")
    sys.stdout.flush()


# ---------------------------------------------------------------------------
# Engine dispatch
# ---------------------------------------------------------------------------

def _search_hmmer(cut_fa, aa_fa, db_path, n_workers):
    """Six-frame translate windows and exhaustively search. Returns
    (hits, win_lengths, aa_fa) — win_lengths maps window_id -> nucl length."""
    log.info("  Six-frame translating windows")
    if os.path.exists(aa_fa + ".fxi"):
        os.remove(aa_fa + ".fxi")
    win_lengths = translate_fasta(cut_fa, aa_fa)
    hmms = load_hmms(db_path)
    optimized = build_optimized_profiles(hmms, alphabet=AMINO_ALPHABET)
    aa_block = build_sequence_block(aa_fa, AMINO_ALPHABET)
    log.info("  Searching {} models over {} frames".format(
        len(hmms), len(aa_block)))
    hits = legacy_search(hmms, aa_block, optimized=optimized)
    log.info("  {} raw domain hits".format(len(hits)))
    return hits, win_lengths


def _search_bath(cut_fa, db_path, db_name, n_workers, outdir):
    """bathsearch --fs directly on the nucleotide windows. win_lengths is
    unused (BATH tblout coords are already nucleotide)."""
    hits = bath_search.run_and_parse(
        db_path, cut_fa, db_name, n_workers=n_workers, outdir=outdir)
    return hits, {}


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_genome(input_fasta, db_paths, db_alphabets, outdir, prefix,
               use_bath=False, n_workers=4,
               win_size=int(1e6), win_ovl=int(1e5)):
    """Run genome mode end-to-end. db_paths/db_alphabets keyed by db name; only
    amino-acid databases are searched (genome mode finds protein domains)."""
    t_start = time.time()
    engine = "bath" if use_bath else "hmmer"

    aa_dbs = [n for n, a in db_alphabets.items() if a == AMINO_ALPHABET]
    if not aa_dbs:
        raise SystemExit(
            "Genome mode needs at least one amino-acid database (it finds TE "
            "protein domains). DNA-only databases (e.g. sine) are unsupported.")
    skipped = [n for n in db_paths if n not in aa_dbs]
    if skipped:
        log.warning("  Genome mode skips DNA database(s): {}".format(
            ", ".join(skipped)))

    win_size, win_ovl = int(win_size), int(win_ovl)
    if engine == "hmmer" and win_size + win_ovl > _HMMER_MAX_WIN_NUCL:
        new_ovl = min(win_ovl, _HMMER_MAX_WIN_OVL)
        new_size = _HMMER_MAX_WIN_NUCL - new_ovl
        log.warning(
            "  HMMER engine: capping window {}+{} -> {}+{} nt to stay under "
            "pyhmmer's {}-residue limit (overlap still exceeds any TE domain, "
            "so no domains are lost). Use --bath for large windows.".format(
                win_size, win_ovl, new_size, new_ovl, _PYHMMER_MAX_AA))
        win_size, win_ovl = new_size, new_ovl

    cut_fa = os.path.join(outdir, "cut.fa")
    log.info("Windowing genome (win_size={}, win_ovl={}) -> {}".format(
        win_size, win_ovl, cut_fa))
    cut_lengths = cut_genome(input_fasta, cut_fa, win_size=win_size,
                             win_ovl=win_ovl)
    log.info("  {} windows".format(len(cut_lengths)))

    aa_fa = os.path.join(outdir, "{}.windows.aa".format(prefix))

    all_features = []
    for name in aa_dbs:
        db_path = db_paths[name]
        config = DB_CONFIGS.get(name)
        if config is None:
            log.warning("  No classifier config for {}, skipping".format(name))
            continue
        if config["clade_parser"] == "gydb":
            config = dict(config)
            config["_clade_map"] = load_gydb_clade_map()

        log.info("--- Genome search: {} ({}) ---".format(
            name, os.path.basename(db_path)))
        t0 = time.time()
        if engine == "bath":
            hits, win_lengths = _search_bath(
                cut_fa, db_path, name, n_workers, outdir)
        else:
            hits, win_lengths = _search_hmmer(cut_fa, aa_fa, db_path, n_workers)
        log.info("  Search done in {:.1f}s".format(time.time() - t0))

        n_pass = 0
        for h in hits:
            result = _hit_to_feature(h, engine, config, win_lengths)
            if result is not None:
                all_features.append(result[0])
                n_pass += 1
        log.info("  {} domain features pass filters".format(n_pass))

    log.info("Resolving overlapping features across {} raw features".format(
        len(all_features)))
    resolved = group_resolve_overlaps(all_features)
    log.info("  {} features after overlap resolution".format(len(resolved)))

    # Sort for stable output: chromosome, then start.
    resolved.sort(key=lambda x: (x[0], x[3]))

    seq_path = cut_fa if engine == "bath" else aa_fa
    gff, seqf, summ = write_outputs(resolved, engine, outdir, prefix, seq_path)
    log.info("  Domain GFF3:  {}".format(gff))
    log.info("  Domain seqs:  {}".format(seqf))
    log.info("  Summary:      {}".format(summ))
    log.info("Genome mode done in {:.1f}s".format(time.time() - t_start))
    return gff, seqf, summ
