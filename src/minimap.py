"""
minimap.py — minimap2 wrapper for pass-2 similarity search.

Unlike mmseqs `--search-type 3`, minimap2 emits multiple chained alignments per
(query, target) pair, so we can legitimately union query- and target-side
intervals across chains to compute both qcov and tcov. No 500 bp chain-gap
heuristic needed: minimap2's own chaining already groups near-diagonal
minimizer seeds, and each PAF row is one such chain. Across-chain union happens
here, in Python.

Filter semantics enforced by the consumer (blast_pass2.classify_from_blast):
  identity ≥ I, qcov ≥ C, tcov ≥ C, length ≥ L   (from --pass2-rule I-C-L)

where qcov = |union of aligned query intervals| / qlen and symmetrically for
tcov.
"""

import logging
import os
import shutil
import subprocess
from collections import OrderedDict, defaultdict

log = logging.getLogger(__name__)


def check_minimap2(bin="minimap2"):
    if shutil.which(bin) is None:
        raise RuntimeError(
            f"{bin!r} not found on PATH. Install minimap2 "
            f"(conda: `mamba install -c bioconda minimap2`) and retry."
        )


def minimap2_version(bin="minimap2"):
    if shutil.which(bin) is None:
        log.warning(f"{bin!r} not found on PATH")
        return None
    r = subprocess.run([bin, "--version"], capture_output=True, text=True, check=False)
    v = (r.stdout or r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr) else "unknown"
    log.info(f"minimap2 version: {v}")
    return v


def run_minimap2(query_fa, target_fa, paf_out, ncpu=4,
                 preset="asm20", extra="",
                 minimap2_bin="minimap2"):
    """
    Run minimap2 and write PAF (with CIGAR / AS tags via -c).

    Defaults:
      -x {preset}    e.g. asm20 for ~20% divergence LTR-RTs
      -c             emit CIGAR + AS:i: alignment score tag (needed for ranking)
      -N 50          keep up to 50 secondary alignments per query
      -p 0.1         accept secondaries down to 10% of primary score
      --secondary=yes
      -I 100G        don't split index (fine for LTR libraries up to a few GB)
      -t ncpu
    """
    os.makedirs(os.path.dirname(os.path.abspath(paf_out)) or ".", exist_ok=True)

    cmd = (
        f"{minimap2_bin} -c -x {preset} -N 50 -p 0.1 --secondary=yes "
        f"-I 100G -t {ncpu} {extra} "
        f"-o {paf_out} {target_fa} {query_fa}"
    )
    log.info(f"minimap2 cmd: {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    if r.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed (exit {r.returncode}): "
            f"{(r.stderr or r.stdout)[:2000]}"
        )
    return paf_out


class PAFRecord:
    """One minimap2 PAF line. 0-based half-open coords on query and target."""
    __slots__ = (
        "qseqid", "qlen", "qstart", "qend", "strand",
        "sseqid", "tlen", "tstart", "tend",
        "matches", "alnlen", "mapq", "score",
    )

    def __init__(self, line):
        vals = line.rstrip("\n").split("\t")
        self.qseqid = vals[0]
        self.qlen = int(vals[1])
        self.qstart = int(vals[2])
        self.qend = int(vals[3])
        self.strand = vals[4]
        self.sseqid = vals[5]
        self.tlen = int(vals[6])
        self.tstart = int(vals[7])
        self.tend = int(vals[8])
        self.matches = int(vals[9])
        self.alnlen = int(vals[10])
        self.mapq = int(vals[11])
        # AS tag = Smith-Waterman alignment score (when -c is used)
        self.score = None
        for tag in vals[12:]:
            if tag.startswith("AS:i:"):
                self.score = int(tag[5:])
                break
        if self.score is None:
            # fallback: use matches as score proxy
            self.score = self.matches


def _union_len(intervals):
    """Length of the union of [lo, hi) half-open intervals."""
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    total += cur_hi - cur_lo
    return total


def parse_paf_besthit(paf_path):
    """Parse PAF, union per (query, target) pair on both axes, return one record
    per (qseqid, sseqid) pair with merged qcov/tcov/identity; then best-hit-per-
    query is selected by caller.

    Returns: list of merged hit dicts with keys:
      qseqid, sseqid, qlen, tlen,
      qcov (0..1), tcov (0..1), fident (0..1),
      alnlen (union on query axis, i.e. unique query bases aligned),
      matches, score  (max AS across chains for this pair)
    """
    if not os.path.exists(paf_path) or os.path.getsize(paf_path) == 0:
        return []

    groups = defaultdict(list)
    with open(paf_path) as f:
        for line in f:
            if not line.strip():
                continue
            r = PAFRecord(line)
            groups[(r.qseqid, r.sseqid)].append(r)

    merged = []
    for (q, t), recs in groups.items():
        q_intervals = [(r.qstart, r.qend) for r in recs]
        t_intervals = [(r.tstart, r.tend) for r in recs]
        q_union = _union_len(q_intervals)
        t_union = _union_len(t_intervals)

        qlen = recs[0].qlen
        tlen = recs[0].tlen

        total_aln = sum(r.alnlen for r in recs)
        total_matches = sum(r.matches for r in recs)
        fident = (total_matches / total_aln) if total_aln else 0.0

        best_score = max(r.score for r in recs)

        merged.append({
            "qseqid": q,
            "sseqid": t,
            "qlen": qlen,
            "tlen": tlen,
            "qcov": q_union / qlen if qlen else 0.0,
            "tcov": t_union / tlen if tlen else 0.0,
            "fident": fident,
            "alnlen": q_union,
            "matches": total_matches,
            "score": best_score,
        })
    return merged


def besthit_per_query(merged):
    """Return OrderedDict[qseqid -> best merged record by score]."""
    best = OrderedDict()
    # Stable sort by score desc for determinism.
    merged_sorted = sorted(merged, key=lambda m: -m["score"])
    for m in merged_sorted:
        if m["qseqid"] not in best:
            best[m["qseqid"]] = m
    return best
