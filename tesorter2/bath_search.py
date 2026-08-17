"""
bath_search.py — Run the BATH aligner as a drop-in replacement for the
pyhmmer/HMMER search on amino-acid databases.

BATH (bathsearch) performs a frameshift-aware translated search of protein
HMMs directly against *nucleotide* sequences, so unlike the HMMER path it does
NOT need six-frame translation of the input. This module:

  1. Resolves a BATH-format HMM for a given repo HMMER3 database (reusing a
     cached / pre-converted .bath.hmm when available, else running bathconvert).
  2. Runs bathsearch (--fs) against the raw nucleotide FASTA.
  3. Parses the tblout into hit dicts matching the schema that
     results.store_legacy expects, so the existing classifier runs unchanged.

BATH's tblout has no posterior-probability ("acc") column, so acc is set to
1.0 — the min_acc>=0.5 filter becomes a no-op for BATH hits, and effective
filtering happens on E-value / coverage / normalized score (the same fields
TEsorter's plant benchmark filtered on).
"""

import logging
import os
import shutil
import subprocess
import tempfile

log = logging.getLogger(__name__)

_RELEASES_URL = "https://github.com/TravisWheelerLab/BATH/releases"

# Fallback for the original development environment. Only consulted when
# BATH_BIN_DIR is unset and the binary is not on PATH.
_LEGACY_BATH_BIN_DIR = "/anvil/projects/x-bio250374/daniel/software/BATH_2/BATH/src"


def _bin(name):
    """Absolute path to a BATH binary.

    Resolution order: BATH_BIN_DIR (an explicit setting always wins, and fails
    loudly rather than silently falling through), then PATH, then the legacy
    development directory.
    """
    bin_dir = os.environ.get("BATH_BIN_DIR")
    if bin_dir:
        path = os.path.join(bin_dir, name)
        if os.path.isfile(path):
            return path
        raise FileNotFoundError(
            f"BATH binary '{name}' not found in BATH_BIN_DIR ({bin_dir}).")

    found = shutil.which(name)
    if found:
        return found

    path = os.path.join(_LEGACY_BATH_BIN_DIR, name)
    if os.path.isfile(path):
        return path

    raise FileNotFoundError(
        f"BATH binary '{name}' not found on PATH. BATH is not packaged on "
        f"conda: download a prebuilt binary from {_RELEASES_URL} (or build "
        f"from source), then put it on PATH or set BATH_BIN_DIR.")


def require_binaries():
    """Fail early if BATH is missing.

    Genome mode searches protein profiles exclusively with BATH, so this runs
    before any expensive setup rather than letting the run window a genome and
    then die at the first search. Re-raised as SystemExit: a missing external
    dependency is a user-facing configuration problem, not a crash, so it
    should print the instructions rather than a traceback.
    """
    for name in ("bathsearch", "bathconvert"):
        try:
            _bin(name)
        except FileNotFoundError as exc:
            raise SystemExit(str(exc))


# Optional pre-converted BATH HMMs, keyed by TEsorter2 db alias. Used as a
# shortcut to skip bathconvert when these files already exist; the cache /
# fallback conversion below is the self-contained path that always works.
_PROJECT = "/anvil/projects/x-bio250374/daniel/bath_vs_hmmer"
KNOWN_CONVERTED = {
    "rexdb": os.path.join(_PROJECT, "analysis_for_paper/hmm_combine/hmm_dbs/rexdb.bath.hmm"),
    "gydb":  os.path.join(_PROJECT, "analysis_for_paper/hmm_combine/hmm_dbs/gydb.bath.hmm"),
    "line":  os.path.join(_PROJECT, "analysis_for_paper/bath_dbs/line_kapitonov.bath.hmm"),
    "tir":   os.path.join(_PROJECT, "analysis_for_paper/bath_dbs/tir_pnas.bath.hmm"),
}


def resolve_bath_db(hmm_path, db_name=None):
    """Return a path to a BATH-format HMM for the given repo HMMER3 database.

    Resolution order:
      1. Cached sidecar {hmm_path}.bath.hmm (same pattern the repo uses for
         .similarity_graph.db) — reuse if present.
      2. A known pre-converted file for this alias (KNOWN_CONVERTED).
      3. Run bathconvert on hmm_path, writing the cache, and use that.
    """
    cache = hmm_path + ".bath.hmm"
    if os.path.isfile(cache):
        log.info(f"  BATH HMM (cached): {cache}")
        return cache

    if db_name and db_name in KNOWN_CONVERTED and os.path.isfile(KNOWN_CONVERTED[db_name]):
        log.info(f"  BATH HMM (pre-converted): {KNOWN_CONVERTED[db_name]}")
        return KNOWN_CONVERTED[db_name]

    log.info(f"  Converting {os.path.basename(hmm_path)} to BATH format -> {cache}")
    # Write to a temp file first so an interrupted convert never leaves a
    # half-written cache that a later run would happily reuse.
    fd, tmp = tempfile.mkstemp(suffix=".bath.hmm", dir=os.path.dirname(cache) or ".")
    os.close(fd)
    try:
        subprocess.run([_bin("bathconvert"), tmp, hmm_path], check=True)
        shutil.move(tmp, cache)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    return cache


def run_bathsearch(bath_hmm, nucl_fasta, tblout_path, n_workers=4,
                   frameshift=True, report_evalue=0.01):
    """Run bathsearch, writing a parseable tblout. Returns tblout_path.

    frameshift=True enables --fs (BATH's frameshift-aware mode), which is the
    whole reason to prefer BATH over HMMER for degraded elements.
    report_evalue caps the reported hits to keep the tblout small; the real
    TEsorter cutoff (1e-3) is applied later by the classifier's filters.
    """
    cmd = [_bin("bathsearch")]
    if frameshift:
        cmd.append("--fs")
    cmd += ["--cpu", str(max(1, n_workers)),
            "-E", str(report_evalue),
            "--tblout", tblout_path,
            "-o", os.devnull,
            bath_hmm, nucl_fasta]
    log.info(f"  bathsearch: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return tblout_path


def parse_bath_tblout(tblout_path):
    """Parse a BATH tblout into hit dicts for results.store_legacy.

    Column layout, relative to `base` (the index of the target-name column):
      +0 target  +2 query  +4 hmm_len  +5 hmm_from  +6 hmm_to  +7 seq_len
      +8 ali_from  +9 ali_to  +10 env_from  +11 env_to  +12 E-value  +13 score
      +14 bias  ...

    BATH 2.0 prepends a numeric "hit ID" column (base=1) and inserts a PID
    column after bias; BATH 1.x has neither (base=0). We never read past +14
    (bias), so the differing trailing columns are irrelevant -- only `base`
    matters. `base` is read once from the column-name header line, which is
    always emitted: it starts with "hit ID" for 2.0 and "target name" for 1.x.
    (Detecting `base` from the data rows would be ambiguous: a purely-numeric
    target name -- e.g. an Ensembl chromosome "1" -- is indistinguishable from
    a 2.0 numeric hit-ID.)
    """
    hits = []
    base = 0
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                header = line.lstrip("#").split()
                if header[:2] == ["hit", "ID"]:
                    base = 1
                elif header[:2] == ["target", "name"]:
                    base = 0
                continue
            p = line.split()
            if len(p) < base + 15:
                continue

            seq_name = p[base + 0]
            hmm_name = p[base + 2]
            hmm_len = int(p[base + 4])
            hmm_from = int(p[base + 5])
            hmm_to = int(p[base + 6])
            seq_len = int(p[base + 7])
            ali_from = int(p[base + 8])
            ali_to = int(p[base + 9])
            env_from = int(p[base + 10])
            env_to = int(p[base + 11])
            evalue = float(p[base + 12])
            score = float(p[base + 13])
            bias = float(p[base + 14])

            # BATH reports descending coords on the minus strand. Encode strand
            # in the target suffix (matching results._parse_frame_info's
            # |fwd1/|rev1 convention) and normalize coords to ascending so the
            # downstream overlap / ordering geometry matches the HMMER path.
            strand_suffix = "fwd1" if ali_from <= ali_to else "rev1"
            ali_lo, ali_hi = sorted((ali_from, ali_to))
            env_lo, env_hi = sorted((env_from, env_to))

            hits.append({
                "target_name": f"{seq_name}|{strand_suffix}",
                "target_len": seq_len,
                "query_name": hmm_name,
                "query_len": hmm_len,
                "evalue": evalue,
                "score": score,
                "bias": bias,
                "dom_num": 1,
                "dom_of": 1,
                "c_evalue": evalue,
                "i_evalue": evalue,
                "dom_score": score,
                "dom_bias": bias,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_lo,
                "ali_to": ali_hi,
                "env_from": env_lo,
                "env_to": env_hi,
                # BATH tblout has no posterior column; 1.0 makes the acc filter
                # a no-op (see module docstring).
                "acc": 1.0,
            })
    return hits


def run_and_parse(hmm_path, nucl_fasta, db_name, n_workers=4, outdir=None,
                  frameshift=True):
    """Resolve the BATH HMM, run bathsearch, and return parsed hit dicts."""
    bath_hmm = resolve_bath_db(hmm_path, db_name=db_name)

    tmp_tblout = None
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        tblout = os.path.join(outdir, f"{db_name}.bath.tblout")
    else:
        fd, tblout = tempfile.mkstemp(suffix=".bath.tblout")
        os.close(fd)
        tmp_tblout = tblout

    try:
        run_bathsearch(bath_hmm, nucl_fasta, tblout, n_workers=n_workers,
                       frameshift=frameshift)
        hits = parse_bath_tblout(tblout)
    finally:
        if tmp_tblout and os.path.exists(tmp_tblout):
            os.remove(tmp_tblout)
    log.info(f"  BATH: {len(hits)} hits for {db_name}")
    return hits
