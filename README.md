# TEBinSorter

Near-perfect replication of [TEsorter](https://github.com/zhangrengang/TEsorter) at >500x speed. Currently an incomplete replication of the original: this code stops at the end of the HMMsearch phase it seeks to accelerate.

## How it works

TEsorter classifies transposable elements by searching translated sequences against HMM profile databases using HMMER's `hmmscan` with `--nobias`. This disables fast filtering for every sequence-model comparison, making the vast majority of runtime a waste: obvious non-hits are screened at full cost so that a handful of true positives aren't missed.

TEBinSorter replaces this with in-process `hmmsearch` via [pyhmmer](https://pyhmmer.readthedocs.io/), IQR-balanced parallel workload scheduling, and an optional **facet classification** mode that uses spliced sub-HMMs for fast pre-screening.

Results are stored in a SQLite database, enabling post-hoc filtering and analysis without re-running searches.

## Search modes

### Default mode (exact TEsorter replication)

Single-pass nobias `hmmsearch` against all models. Produces identical results to TEsorter but dramatically faster via:
- `hmmsearch` instead of `hmmscan` (no per-sequence database reload)
- In-process pyhmmer instead of subprocess calls
- IQR-based parallel strategy: models with outlier M² cost use `parallel='targets'`, the rest use `parallel='queries'`, avoiding thread starvation from large models

### Facet mode (`--facet`)

For amino acid databases, uses spliced sub-HMMs ("facets") for fast pre-screening:

1. **Facet screen**: Tiered sub-HMMs (96→64→48→32 positions) searched against all frames. Uniformly sized, perfect parallel load balance.
2. **Targeted verification**: Top facet hit per domain family verified with single full-model nobias search.
3. **Cross-family completion**: Classified frames searched for missing domain families.
4. **Legacy fallback**: Unclassified frames (no facet signal) get full nobias search.

DNA databases (AnnoSINE) always use the default legacy search -- DNA facets do not provide sufficient sensitivity gains to justify the overhead.

## Benchmarks

Rice TE library (2,431 sequences), 4 processors.

**Important caveat:** TEsorter baseline times were measured on WSL2 (Windows Subsystem for Linux), which has known I/O and process-spawning overhead compared to native Linux. TEsorter's architecture (many subprocess calls to hmmscan, heavy temp file I/O) is disproportionately affected by this. The true speedup on native Linux may be smaller. TEBinSorter times are also WSL2 but its architecture (in-process pyhmmer, minimal disk I/O) is less sensitive to this penalty.

### Default mode

| Database | Models | TEsorter | TEBinSorter | Speedup |
|----------|--------|----------|-------------|---------|
| TIR | 17 | 26m 41s | 0.9s | ~1,779x |
| LINE | 28 | 56m 45s | 1.7s | ~2,003x |
| SINE | 88 | ~30m | 38.1s | ~47x |
| REXdb | 266 | >2h | 20.2s | >356x |
| GyDB | 314 | >2h | 18.3s | >393x |

All 5 databases searched in **79 seconds** total (4 processors).

### Facet mode (amino acid databases)

Rice TE library, REXdb (266 models), with cross-family completion:

| Metric | Value |
|--------|-------|
| Post-filter family recall vs default | 99.8% |
| Misses | 4 out of 2,090 |
| Facet build + screen | 11.4s |
| Verification | 2.4s |
| Cross-family completion | 3.3s |
| Legacy fallback | 1.2s |

### Verification against TEsorter

Tested on the rice6.9.5.liban TE library (2,431 sequences). Default mode results filtered using TEsorter's exact thresholds (`coverage >= 20`, `evalue <= 1e-3`, `probability >= 0.5`, `normalized score >= 0.1`) and compared at the (sequence, model) pair level.

| Database | TEsorter pairs | TEBinSorter pairs | Common | TEsorter only | Notes |
|----------|---------------|-------------------|--------|---------------|-------|
| TIR | 430 | 430 | 430 | 0 | Perfect match |
| LINE | 806 | 778 | 778 | 28 | TEsorter rounding artifact |
| REXdb | 2,787 | 2,784 | 2,767 | 20 | 17 rounding, 3 bias filter |
| GyDB | 7,496 | 7,448 | 7,448 | 48 | All rounding artifacts |
| SINE | 0 (TEsorter bug) | 731 | -- | -- | TEsorter forces `seq_type='prot'` for a DNA database |

**Alignment coordinates:** Identical to TEsorter's GFF3 output (verified 430/430 on TIR).

**E-values:** Identical (Z set to match hmmscan convention).

**TEsorter rounding bug:** TEsorter rounds `domain_score / model_length` to 2 decimal places before threshold comparison. TEBinSorter uses full precision by default. `--compat-tesorter-rounding` replicates the old behavior.

## Installation

### Dependencies

- Python >= 3.10
- [pyhmmer](https://pyhmmer.readthedocs.io/) >= 0.10
- [pyfastx](https://github.com/lmdu/pyfastx) >= 2.0
- numpy

```bash
pip install pyhmmer pyfastx numpy
```

### Clone

```bash
git clone https://github.com/KGerhardt/TEBinSorter.git
cd TEBinSorter
```

## Usage

### Basic search (default mode)

```bash
python3 src/pipeline.py input.fasta -d rexdb -p 4 -o output_dir
```

### Multiple databases

```bash
python3 src/pipeline.py input.fasta -d rexdb,gydb,tir -p 4 -o output_dir
```

### All databases

```bash
python3 src/pipeline.py input.fasta --max-search -p 4 -o output_dir
```

### Facet mode (faster, AA databases only)

```bash
python3 src/pipeline.py input.fasta -d rexdb --facet -p 4 -o output_dir
```

DNA databases automatically fall back to default mode when `--facet` is specified.

### Emit partitions for BATH aligner

```bash
python3 src/pipeline.py input.fasta -d rexdb --emit-bath -o output_dir
# Partitions written to output_dir/BATHwater/
```

Or standalone from an existing results database:
```bash
python3 src/emit.py results.db sequences.aa database.hmm -o output_dir/BATHwater
```

## Databases

TEBinSorter auto-detects database alphabet (amino acid or DNA) from the HMM file. Built-in aliases:

| Alias | Database | Alphabet |
|-------|----------|----------|
| `rexdb` | REXdb v4.0 + metazoa v3.1 | amino |
| `gydb` | GyDB 2.0 | amino |
| `line` | Kapitonov et al. LINE RT | amino |
| `tir` | Yuan & Wessler TIR TPase | amino |
| `sine` | AnnoSINE | DNA |

Custom HMM databases can be passed as file paths in the `-d` argument. Pre-computed similarity graphs are generated automatically on first use and cached alongside the HMM file.

## Output files

| File | Description |
|------|-------------|
| `{prefix}.db` | SQLite database with all results |
| `{prefix}.aa` | Six-frame translated amino acid sequences (indexed) |
| `{prefix}.{db}.classifications.tsv` | Facet classifications with confidence tiers (facet mode) |
| `BATHwater/` | Routed FASTA partitions (with `--emit-bath`) |

## Architecture

### Core modules

- **`pipeline.py`** — Main CLI and search orchestration
- **`search.py`** — HMM search engine with IQR-balanced parallelism
- **`sequence.py`** — FASTA ingestion (pyfastx) and six-frame translation (pyhmmer)
- **`hmm.py`** — HMM loading, alphabet detection, optimized profile construction
- **`results.py`** — SQLite persistence with pre-parsed columns (base_seq, strand, frame, domain_type)
- **`deconflict.py`** — Numpy-based hit deconfliction and parameterized filtering

### Facet classification

- **`decompose_hmm.py`** — Standalone sub-HMM decomposition and splicing. Tiered window sizes, configurable overlap. Works with DNA and amino acid HMMs.
- **`model_graph.py`** — Cross-model similarity graph. Pre-built for all databases.
- **`facet_classify.py`** — Facet screen → verification → cross-family completion → legacy fallback
- **`cross_family.py`** — Targeted search for missing domain families in classified frames

### Utilities

- **`emit.py`** — BATH aligner partition emission (standalone or pipeline-integrated)
- **`id_registry.py`** — Deterministic integer ID registry for fast numeric storage
