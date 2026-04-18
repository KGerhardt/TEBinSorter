# TEBinSorter

Near-perfect replication of [TEsorter](https://github.com/zhangrengang/TEsorter) at greatly improved speed. Currently an incomplete replication of the original: this code stops at the end of the HMMsearch phase it seeks to accelerate.

## How it works

TEsorter classifies transposable elements by searching translated sequences against HMM profile databases using HMMER's `hmmscan` with `--nobias`. This disables fast filtering for every sequence-model comparison, making the vast majority of runtime a waste: obvious non-hits are screened at full cost so that a handful of true positives aren't missed.

TEBinSorter achieves acceleration through substantial architectural changes to almost every portion of the code.

* Rather than HMMscan, TEBinSorter utilizes `hmmsearch` via [pyhmmer](https://pyhmmer.readthedocs.io/), which is generally more performant. E-value differences between HMMscan and HMMsearch are eliminated by parameterization to ensure identical results, just faster.
* TEBinSorter utilizes intelligent parallel workload balancing to optimize the efficiency of searches and ensure near-perfect CPU utilization
* TEBinSorter implements an alternative, optional algorithm for more rapidly searching protein databases (most of what TESorter uses) by pre-screening TEs to identify probable best hits before expensive searches are performed.

Results are stored in a SQLite database, enabling post-hoc filtering and analysis without re-running searches.

## Search modes

### Default mode (exact TEsorter replication)

Single-pass nobias `hmmsearch` against all models. Produces identical results to TEsorter but dramatically faster via:
- `hmmsearch` instead of `hmmscan` (no per-sequence database reload)
- In-process pyhmmer instead of subprocess calls
- HMM model cost-aware parallel load balancing to prevent idle CPU time

### Facet mode (`--facet`)

For amino acid databases, uses spliced sub-HMMs ("facets") for fast pre-screening:

1. **Facet screen**: Tiered sub-HMMs (96→64→48→32 amino acids) searched against all six translated protein reading frames of each input TE. Uniformly sized, perfect parallel load balance is intrinsic.
2. **Targeted verification**: Top facet hit per domain family verified with single full-model nobias search.
3. **Cross-family completion**: Classified frames searched for missing domain families.
4. **Legacy fallback**: Unclassified frames (no facet signal) optionally get full nobias search.

DNA databases (AnnoSINE) always use the default legacy search -- DNA facets do not provide sufficient sensitivity gains to justify the overhead.

## Quick notes ##

* All results reported by TEBinSorter ultimately emerge from identical HMM alignments to identical sequences with identical locations and hit probability metrics (e.g. E-values) using the same --nobias logic as the original; this completely guarantees that a hit found by TEBinSorter is bitwise identical to a hit found by TESorter.
* In the facets mode, TEBinSorter is finding a relatively small subset of all TESorter hits that are extremely likely to be the single best hit per input sequence. In a sense, this means that the facets search finds less. However, after TESorter filters results, the filtering produces a subset of hits essentially identical (99.98%) to those found by TEBinSorter - TEBinSorter finds these hits directly instead of by post-processing from a completely exhaustive search and skips the cost of the exhaustive search on most sequences.

## Benchmarks

Rice TE library (2,431 sequences), 4 processors.

**Important caveat:** TEsorter baseline times were measured on WSL2 (Windows Subsystem for Linux), which has known I/O and process-spawning overhead compared to native Linux. TEsorter's architecture (many subprocess calls to hmmscan, heavy temp file I/O) is disproportionately affected by this. The true speedup on native Linux may be smaller. TEBinSorter times are also WSL2 but its architecture (in-process pyhmmer, minimal disk I/O) is less sensitive to this penalty. I'm going to measure a true comparsion a bit later - the numbers will not be nearly as impressive as they appear here, unless you too are running the program on WSL.

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
| LINE | 806 | 778 | 778 | 28 | All rounding artifacts |
| REXdb | 2,787 | 2,784 | 2,767 | 20 | 17 rounding, 3 bias filter |
| GyDB | 7,496 | 7,448 | 7,448 | 48 | All rounding artifacts |
| SINE | 0 (TEsorter bug) | 731 | -- | -- | TEsorter forces `seq_type='prot'` for a DNA database |

**Alignment coordinates:** Identical to TEsorter's GFF3 output

**E-values:** Identical (Z set to match hmmscan convention).

**TEsorter rounding bug:** TEsorter rounds `domain_score / model_length` to 2 decimal places before threshold comparison. This results in cases where 9.9 rounds -> 10 and therefore passes a TESorter filter it never should have. TEBinSorter uses full precision by default. `--compat-tesorter-rounding` replicates the old rounding behavior.

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

### Emit partitions for BATH aligner - currently in need of revision.

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

# Extended methods

## HMM facets

The concept of an HMM facet is one I've been exploring for other projects but found a use for in this one. By utlizing pyhmmer, we can observe the emission probabilities of an HMM model in a friendly manner in Python. I use these emission probabilities to find conserved subregions of the HMM model which are highly influential in the decisionmaking process of HMMer to actually search a sequence. These regions are extracted to form a "facet." There are multiple advantages to this approach:

* Short models tend towards specificity. A full-length model will make an effort to compile information about required search effort across an entire sequence, while short models confirm or reject a local region quickly.
* Much of the information that the full-length model would glean about the appropriateness of a a search from a long view of a sequence is extremely highly correlated with what a good facet will report. A good facet hit almost ensures a good full-length model hit.
* Facets can be intelligently sized so that they exactly pack SIMD lanes (96, 64, 32 amino acids, get consumed mostly in 16-AA sized bites) in the HMMer internals. This reduces low-level CPU waste compared to less politely divisible sizes of sequence.

Are they mathematically correct or statistically sound? I honestly have no idea. But they work.
