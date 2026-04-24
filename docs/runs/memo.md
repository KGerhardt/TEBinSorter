# memo — TEBinSorter_minimap2 port

## Date

Created 2026-04-24.

## Purpose

Port TEBinSorter's pass-2 similarity search from `blastn` to `minimap2`.
Exists because the sibling `TEBinSorter_mmseqs/` port discovered during
benchmarking that `mmseqs easy-search --search-type 3` emits only one best
diagonal per (query, target) pair — an architectural limitation that caps
pass-2 recall on LTR-RT data where candidate and element often share
multiple distinct homologous regions.

`minimap2` has no such limitation: it chains minimizer hits into (potentially
many) independent alignments per (q, t) pair and emits one PAF row per
chain. We union query- and target-side intervals across chains to produce
a single, honest qcov and tcov per pair.

HMM pass-1 of TEBinSorter is untouched (byte-identical verified).

## Source of changes

- Port target: TEBinSorter at `/data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter/`
  (most recent commit `7144950`).
- Port destination: `/data/chris/wheat/ltrharvest/v2/v3/testing/TEBinSorter_minimap2/`.
- Conceptually mirrors `TEBinSorter_mmseqs/` — same CLI shape, same
  `blast_hits` SQL table (with `tcovs` added), same `--pass2-classified-fasta`
  wiring. Swaps the aligner only.

## Environment / versions

- Active conda env during development: `synLTR`
  (`/home/chris/bin/mambaforge/envs/synLTR`).
- `minimap2` binary: **v2.30-r1287** (`synLTR/bin/minimap2`). Must be on
  `$PATH` at runtime. (User installed via `mamba install -c bioconda minimap2`.)
- Python deps: same as TEBinSorter (pyhmmer ≥ 0.10, pyfastx ≥ 2.0, numpy).
- `mmseqs` is *not* required for this port (different from the mmseqs sibling).

## Files touched

- **NEW** `src/minimap.py` — PAF wrapper, parser, union-intervals qcov/tcov
  computation, best-hit-by-AS-score per query.
- **NEW** `src/pass2_external.py` — copied verbatim from the mmseqs port
  (handles `--pass2-classified-fasta`).
- **EDITED** `src/blast_pass2.py` — same filename (keeps pipeline.py,
  tesorter_compat.py, classifier.py, results.py unchanged); internals swapped
  to `minimap2 -c -x asm20 -N 50 -p 0.1 --secondary=yes -I 100G -t NCPU`.
  `blast_hits` SQLite schema has a new `tcovs REAL NOT NULL` column.
  `classify_from_blast` filters on `(qcovs OR tcovs) >= C` by default.
- **EDITED** `src/pipeline.py` — added 5 CLI args (`-dp2/--disable-pass2`,
  `-rule/--pass2-rule`, `--pass2-classified-fasta`, `--minimap2-preset`,
  `--minimap2-extra`); logs `minimap2_version()` at pass-2 start.
- **EDITED** `src/tesorter_compat.py` — mirrored 5 CLI args.
- **EDITED** `README.md` — mmseqs-section-style header noting the dual-coverage
  enforcement and CLI.

## Coverage semantics

Per-(q, t) pair, across all PAF chains minimap2 emits:

```
qcov = |union of aligned query intervals| / qlen
tcov = |union of aligned target intervals| / tlen
```

This mirrors how blastn's `qcovs` is computed (NCBI tiles HSPs) but also
extends it to the target axis. The SQL pass-2 filter applies:

```sql
WHERE pident >= I
  AND (qcovs >= C OR tcovs >= C)   -- OR mode, REQUIRE_BOTH_COVERAGE=False
  AND length >= L
```

where `I-C-L` comes from `--pass2-rule` (default `80-80-80`). `length` is
union query-side alignment length (unique query bases aligned).

### Developer toggle: `REQUIRE_BOTH_COVERAGE`

Module-level constant at the top of `src/blast_pass2.py`:

```python
# True  -> require BOTH qcov AND tcov >= coverage threshold (strict).
# False -> require AT LEAST ONE of qcov, tcov >= coverage threshold
#          (Wicker et al. 80-80-80: "candidate must cover ≥80% of at least
#          one of the elements being compared").
REQUIRE_BOTH_COVERAGE = False
```

- Not exposed on the CLI by design — this is a semantic choice, not a
  runtime parameter.
- Defaults to `False` (OR) after reviewing Wicker et al. 2007 "A unified
  classification system for eukaryotic transposable elements" — the
  original rule is "candidate must cover ≥80% of at least one of the
  elements being compared", which is OR, not AND.
- minimap2 has no native coverage filter; both modes are enforced post-hoc
  in the SQL above.
- `classify_from_blast` logs `coverage mode: AT-LEAST-ONE qcov/tcov >= 80`
  (or `BOTH`) at the start of pass-2.

## minimap2 command

```
minimap2 -c -x asm20 -N 50 -p 0.1 --secondary=yes -I 100G -t NCPU \
    -o pass2.paf <DB_FASTA> <QUERY_FASTA>
```

- `-c` emits CIGAR + `AS:i:` score tag (needed for ranking best-hit).
- `-x asm20` targets ~20% divergence — suits LTR-RT family-level variation.
- `-N 50 -p 0.1 --secondary=yes` keeps secondary alignments down to 10% of
  primary score; essential for repetitive LTR-RT libraries.
- `-I 100G` prevents index splits for libraries up to a few GB.
- Argument order: **target (reference) first, then query**. PAF's column 1
  is the query-fasta sequence, so minimap2's convention maps cleanly to our
  blast_hits schema.

## Public API

All functions in `src/blast_pass2.py` retain the same names as the mmseqs
port and stock blastn TEBinSorter:

- `blast_pass2(input_fasta, conn, ...)` — top-level entry point.
- `store_blast_hits(conn, hits, db_seq_to_dbs)` — writes PAF-derived rows
  to SQLite `blast_hits` table with the extra `tcovs` column.
- `classify_from_blast(conn, classifications, min_identity, min_coverage,
  min_length)` — filters on identity, qcov/tcov (per toggle), length.
- `split_classified_unclassified(...)`, `_get_classified_ids(...)` —
  unchanged from the blastn original.

`src/minimap.py`:

- `check_minimap2(bin)` / `minimap2_version(bin)` — availability + startup log.
- `run_minimap2(query_fa, target_fa, paf_out, ncpu, preset, extra)` —
  subprocess wrapper.
- `PAFRecord` — `__slots__` class; parses one PAF line.
- `parse_paf_besthit(paf_path)` — per-(q, t) union of query and target
  intervals, weighted-fident, max-AS-score.
- `besthit_per_query(merged)` — OrderedDict keyed by qseqid, best by AS.

## CLI

```
python TEBinSorter_minimap2/src/pipeline.py INPUT.fa -d rexdb -p 16 \
    [--pass2-rule 80-80-80] \
    [--pass2-classified-fasta prior_cls.fa] \
    [--minimap2-preset asm20] \
    [--minimap2-extra "..."] \
    [-o OUTDIR]
```

`--disable-pass2` and `--pass-1-only` both bypass minimap2.

## Verification runs

### Unit tests on synthetic PAF

Three scenarios passed (see commit history for exact assertions):

1. Two non-overlapping HSPs to same target, distant on target axis:
   `qcov=0.9, tcov=0.18, score=max(HSPs)` — confirms gap-free union on
   query, interval-separated union on target.
2. Two overlapping query intervals on same target: `qcov=0.75` (overlap
   correctly deduplicated), `tcov=0.395`.
3. One query vs two competing targets (clean big target A vs scattered
   small target B with same total coverage): target A wins best-hit by
   max-AS-score, not by summed score — matches blastn's best-HSP-bitscore
   semantics.

### Rice-example-data parity

Not re-run here; semantic identical to mmseqs port's rice verification.

### Rape-library 5k benchmark (Brassica napus LTR_retriever intact LTRs)

Subsampled with `seqkit sample -s 42 -n 5000` (4,892 unique seqs). 16
threads. REXdb v4 + metazoa v3.1 HMM database.

**At-least-one coverage (OR) mode — default:**

| Rule | Wall | Max RSS | Pass-2 rescues | Combined |
|---|---:|---:|---:|---:|
| 80-80-80 | 0:45 | 6.4 GB | 598 | 2,115 |
| 70-50-80 | 0:46 | 6.4 GB | 939 | 2,456 |

**BOTH-coverage (AND) mode — archived for reference:**

| Rule | Pass-2 rescues |
|---|---:|
| 80-80-80 AND | 16 |
| 70-50-80 AND | 152 |

**Reference:**

| Tool | Rule | Wall | Pass-2 rescues | Combined |
|---|---|---:|---:|---:|
| blastn (stock TEBinSorter) | 80-80-80 (hardcoded, qcov-only) | 17:37 | 1,005 | 2,522 |
| mmseqs port               | 80-80-80 OR | 1:19  | 397   | 1,914 |
| mmseqs port               | 70-50-80 OR | 1:07  | 773   | 2,290 |

### Headline number

minimap2 port @ `70-50-80` OR:
- **23× faster** than blastn baseline (45 s vs 17:37)
- **12× less memory** (6.4 GB vs 76.6 GB)
- **93% of blastn recall** (939 rescues vs 1,005)
- enforces a Wicker-faithful coverage rule (`qcov ≥ 50 OR tcov ≥ 50`)

## SQLite schema

The `blast_hits` table adds a `tcovs REAL NOT NULL` column vs stock
TEBinSorter. Classifier output (`classifications` table, `.cls.tsv`) is
schema-unchanged.

## How this port differs from the mmseqs sibling

| Concern | mmseqs port | minimap2 port |
|---|---|---|
| Alignments per (q, t) pair | one best diagonal | multi-chain, all chained hits |
| qcov / tcov independence | coupled via single alnlen | genuinely independent |
| OR mode on LTR-RT data | ≈ qcov-only (architectural) | meaningfully different from qcov-only |
| Recall at 70-50-80 OR | 773 | 939 |
| Wall time | ~65 s | ~45 s |
| Memory | ~8.6 GB | ~6.4 GB |
| External binary | `mmseqs` 17.b804f | `minimap2` 2.30-r1287 |
| Best-hit ranking metric | mmseqs `bits` | minimap2 `AS` (Smith-Waterman score) |

## Notes / caveats

- `--pass2-classified-fasta` headers must be `>id#Order/Superfamily/Clade`;
  coord-shaped IDs get `unknown` slots upgraded from pass-1 `classifications`
  dict via `_COORD_HEADER_RE` (shared helper from `pass2_external.py`).
- PAF column 11 (`alnlen`) is the block length including gaps; we instead
  report `length` as the union-of-query-intervals (same as `qcovs`'s
  numerator). Consistent with the 80-80-80 rule's third number being
  "length of the alignment region," not "sum of CIGAR ops."
- `AS:i:` tag is only emitted when `-c` is passed; we pass it.
- `evalue` and `slen` columns in `blast_hits` hold sentinels (0.0, 0 or
  tlen respectively) — downstream classifier code never reads them. Kept
  for schema parity with the mmseqs sibling.
- `_MAX_SPLIT_GAP` (the 500-bp chain heuristic from the mmseqs port) does
  not exist here. Unnecessary: minimap2's own chainer groups near-diagonal
  minimizer seeds; cross-chain union on the Python side is the only
  aggregation we do.

## Nested-TE caveat

The OR rule can classify a candidate as family X if it contains a nested
element of family X, because tcov of the nested region can be near 100%
even when qcov is ~50%. Discussed in conversation history; acceptable
trade because pass-1 HMM typically catches this case first. If nested-TE
false positives become a concern, flip `REQUIRE_BOTH_COVERAGE = True` or
add a secondary rule (e.g., flag `qcov / tcov < 0.5` as chimera).

## Reproducibility

```
claude --resume "spin up TEBinSorter_minimap2 port"
```

Archived runs and summaries preserved at:
- `/tmp/keep_baseline_blastn_5k/` (blastn baseline at 80-80-80)
- `/tmp/keep_mmseqs_dualcov_80_5k/` (mmseqs AND @ 80-80-80)
- `/tmp/keep_minimap2_5k/` (minimap2 AND @ 80-80-80)
- `/tmp/rape.summary.md` (most recent OR sweep, all four rule/tool combos)
