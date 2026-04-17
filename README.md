# TEBinSorter

Near-perfect replication of [TEsorter](https://github.com/zhangrengang/TEsorter) at >500x speed. Currently an incomplete replication of the original: this code stops at the end of the HMMsearch phase it seeks to accelerate.

## How it works

TEsorter classifies transposable elements by searching translated sequences against HMM profile databases using HMMER's `hmmscan` with `--nobias`. This disables fast filtering for every sequence-model comparison, making the vast majority of runtime a waste: obvious non-hits are screened at full cost so that a handful of true positives aren't missed.

TEBinSorter replaces this with a two-pass approach:

1. **Pass 1 (coarse screen):** `hmmsearch` with HMMER's default bias filter enabled. This rejects >95% of sequence-model comparisons almost instantly, identifying only the plausible pairs.
2. **Pass 2 (sensitive search):** `hmmsearch` with bias filter disabled, but only against the model subset each sequence was routed to in pass 1. Work is distributed across processes using M²-aware cost scheduling to balance large and small HMM models.

Results are stored in a SQLite database with separate tables for each pass, enabling inspection and filtering without re-running searches.

## Benchmarks

Rice TE library (2,431 sequences), 4 processors.

**Important caveat:** TEsorter baseline times were measured on WSL2 (Windows Subsystem for Linux), which has known I/O and process-spawning overhead compared to native Linux. TEsorter's architecture (many subprocess calls to hmmscan, heavy temp file I/O) is disproportionately affected by this. The true speedup on native Linux may be smaller. TEBinSorter times are also WSL2 but its architecture (in-process pyhmmer, minimal disk I/O) is less sensitive to this penalty.

### Two-pass mode (default)

| Database | Models | TEsorter | TEBinSorter | Speedup |
|----------|--------|----------|-------------|---------|
| TIR (Yuan & Wessler) | 17 | 26m 41s | 2.1s | ~762x |
| LINE (Kapitonov) | 28 | 56m 45s | 4.9s | ~695x |
| SINE (AnnoSINE) | 88 | ~30m | 4.3s | ~419x |
| REXdb (v4 + metazoa) | 266 | >2h | 39.9s | >180x |
| GyDB | 314 | >2h | 44.5s | >162x |

Some of these benchmarks are dominated by the model-loading and sequence routing costs associated with the very small test library of TE sequences used; at larger scales, routing costs amortize and performance gains from a reduced search space become much more obvious. This is most visible in the SINE data for this dataset, as this model contains vastly more expensive individual HMM search operations that dominate runtime in spite of the small target set. The reduction of query sequences produces the ~10x better performance relative to the legacy mode. Similar gains are expected for larger datasets on other sequences.

### Legacy mode (`--legacy`, exact TEsorter replication)

| Database | Models | TEsorter | TEBinSorter | Speedup |
|----------|--------|----------|-------------|---------|
| TIR | 17 | 26m 41s | 0.9s | ~1,779x |
| LINE | 28 | 56m 45s | 1.7s | ~2,003x |
| SINE | 88 | ~30m | 38.1s | ~47x |
| REXdb | 266 | >2h | 20.2s | >356x |
| GyDB | 314 | >2h | 18.3s | >393x |

All 5 databases searched in **79 seconds** total (legacy mode, 4 processors).

By adopting the original no-filter algorithm, `--legacy` mode ensures identical results to the original TEsorter. However, TEBinSorter contains other algorithmic and parallel workload balancing improvements that nonetheless substantially improve overall performance.

### Verification against TEsorter

Tested on the rice6.9.5.liban TE library (2,431 sequences). For each database, TEBinSorter pass-2 results were filtered using TEsorter's exact thresholds (`coverage >= 20`, `evalue <= 1e-3`, `probability >= 0.5`, `normalized score >= 0.1`) and compared at the (sequence, model) pair level.

| Database | TEsorter pairs | TEBinSorter pairs | Common | TEsorter only | Notes |
|----------|---------------|-------------------|--------|---------------|-------|
| TIR | 430 | 430 | 430 | 0 | Perfect match |
| LINE | 806 | 778 | 778 | 28 | All 28 present in raw results; fall below threshold due to TEsorter rounding `score/model_len` to 2 decimal places (e.g. 0.0993 rounds to 0.10) |
| REXdb | 2,787 | 2,784 | 2,767 | 20 | 17 rounding artifacts, 3 genuine bias filter rejections (see below) |
| GyDB | 7,496 | 7,448 | 7,448 | 48 | All 48 are rounding artifacts |
| SINE | 0 (TEsorter bug) | 731 | -- | -- | TEsorter forces `seq_type='prot'` for SINE, which is a DNA database |

**Alignment coordinates:** For all matching pairs, nucleotide coordinates after AA-to-nt conversion are identical to TEsorter's GFF3 output. Verified 430/430 on TIR.

**E-values:** Identical between TEsorter and TEBinSorter for all compared hits (Z set to match hmmscan convention).

**TEsorter rounding bug:** TEsorter rounds `domain_score / model_length` to 2 decimal places *before* comparing against the score threshold. This causes borderline scores like 0.0993 to round up to 0.10 and pass. TEBinSorter uses full precision by default. The `--compat-tesorter-rounding` flag replicates the old behavior for exact reproduction.

**Bias filter sensitivity (3 REXdb misses):** Three sequence-model pairs produce zero signal with the HMMER bias filter enabled but valid hits with it disabled. These sequences have ~5% stop codons in their translated frames -- the bias correction over the full query length masks the genuine domain signal. Setting `--F1 0.1` recovers 2 of 3; the third requires `F1=1.0` (effectively disabling the filter). This is a fundamental limitation of the bias filter approach on heavily non-coding translated sequence, and represents 0.1% of total hits.

## Installation

### Dependencies

- Python >= 3.10
- [pyhmmer](https://pyhmmer.readthedocs.io/) >= 0.10
- [pyfastx](https://github.com/lmdu/pyfastx) >= 2.0
- BLAST+ (planned, not yet required)

```bash
pip install pyhmmer pyfastx
```

### Clone

```bash
git clone https://github.com/yourusername/TEBinSorter.git
cd TEBinSorter
```

## Usage

### Basic search

```bash
python3 src/pipeline.py input.fasta -d rexdb -p 4 -o results
```

### Multiple databases

```bash
python3 src/pipeline.py input.fasta -d rexdb,gydb,tir -p 4 -o results
```

### All databases

```bash
python3 src/pipeline.py input.fasta --max-search -p 4 -o results
```

### Pass 1 only (inspect before committing to pass 2)

```bash
python3 src/pipeline.py input.fasta -d rexdb --pass-1-only -o results
# Inspect results.db with sqlite3 or any SQLite browser
# Then re-run without --pass-1-only when ready
```

### Emit partitions for BATH aligner

From the pipeline:
```bash
python3 src/pipeline.py input.fasta -d rexdb --emit-bath -o results
# Partitions written to results.BATHwater/
```

Or standalone from an existing results database:
```bash
python3 src/emit.py results.db results.aa database.hmm -o results.BATHwater
```

### Tuning the MSV filter

The `--F1` parameter controls the MSV filter stringency in pass 1. Higher values are more permissive, capturing hits in compositionally biased sequences at modest runtime cost:

```bash
python3 src/pipeline.py input.fasta -d rexdb --F1 0.1 -o results
```

Default is 0.02 (HMMER standard). Setting to 0.1 adds ~10% runtime and captures ~3 additional hits per thousand that the default filter rejects.

## Databases

TEBinSorter auto-detects database alphabet (amino acid or DNA) from the HMM file. Built-in aliases:

| Alias | Database | Alphabet |
|-------|----------|----------|
| `rexdb` | REXdb v4.0 + metazoa v3.1 | amino |
| `gydb` | GyDB 2.0 | amino |
| `line` | Kapitonov et al. LINE RT | amino |
| `tir` | Yuan & Wessler TIR TPase | amino |
| `sine` | AnnoSINE | DNA |

Custom HMM databases can be passed as file paths in the `-d` argument.

## Output files

| File | Description |
|------|-------------|
| `{prefix}.db` | SQLite database with all results |
| `{prefix}.aa` | Six-frame translated amino acid sequences (indexed) |
| `{prefix}.pass1.tsv` | Raw pass-1 domain hits |
| `{prefix}.pass2.tsv` | Raw pass-2 domain hits |
| `{prefix}.best.tsv` | Best hit per (sequence, model) pair with nucleotide coordinates |
| `{prefix}.domains.tsv` | All domain occurrences with nucleotide coordinates |
| `{prefix}.domains.faa` | Domain protein subsequences |
| `{prefix}.BATHwater/` | Routed FASTA partitions (with `--emit-bath`) |

## Under construction: facet classification

An experimental faster search mode using **sub-HMM decomposition** ("facets") is in active development. Rather than searching every sequence against every full model, facets are small (64-position) spliced sub-HMMs built from the highest-scoring regions of each model. The approach:

1. **Facet screen**: Search all facets against all frames in one pass. Facets are uniformly small, enabling perfect parallel load balancing.
2. **Targeted verification**: For each frame, verify only the top-scoring facet per domain family with a single full-model nobias search. Typically ~70 unique models instead of 266.
3. **Confidence tiers**: Unverified facet hits are reported with empirically-calibrated confidence levels based on facet score:

| Facet score | Confidence | Empirical verification rate |
|------------|------------|---------------------------|
| >= 20 | confirmed | 100% |
| 15-20 | probable | 99.9% |
| 10-15 | candidate | 86.5% |
| < 10 | excluded | unreliable |

4. **Legacy fallback**: Frames with no facet signal get full nobias legacy search.

### Preliminary results (10k frame sample, REXdb 266 models)

| Metric | Value |
|--------|-------|
| Family-level recall vs legacy | 94.3% |
| Exact best-model match vs legacy | 88.3% |
| Facet screen time | 4.0s |
| Verification time (70 models) | 1.7s |
| Legacy fallback time | 1.4s |
| **Total** | **~13s search time** |

Recall by confidence tier:

| Tier | Assignments | Exact model accuracy |
|------|-------------|---------------------|
| confirmed (>= 20) | 1,760 | 94.7% |
| probable (15-20) | 207 | 85.5% |
| candidate (10-15) | 191 | 67.0% |

### Supporting modules

- **`decompose_hmm.py`** — Standalone sub-HMM decomposition and splicing. Works with DNA and amino acid HMMs. Greedy window selection with configurable overlap, coverage bounds, and score thresholds.
- **`model_graph.py`** — Cross-model similarity graph built from consensus-vs-consensus search. Persisted as SQLite alongside HMM files. Pre-built for all 5 databases.
- **`facet_classify.py`** — Facet screen → targeted verification → legacy fallback pipeline.
- **`iterative_search.py`** — Graph-guided iterative confirmation (alternative to facet classification, higher recall at higher cost).
