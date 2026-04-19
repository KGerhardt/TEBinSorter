# TEBinSorter

Near-perfect replication of [TEsorter](https://github.com/zhangrengang/TEsorter) at greatly improved speed. 

Most differences in TESorter replication are the result of small bugfixes in TEBinSorter.

## How it works

TEsorter classifies transposable elements by searching translated sequences against HMM profile databases using HMMER's `hmmscan` with `--nobias`. This disables fast filtering for every sequence-model comparison, making the vast majority of runtime a waste: obvious non-hits are screened at full cost so that a handful of true positives aren't missed.

TEBinSorter achieves acceleration through substantial architectural changes to almost every portion of the code.

* Rather than HMMscan, TEBinSorter utilizes `hmmsearch` via [pyhmmer](https://pyhmmer.readthedocs.io/), which is generally more performant. E-value differences between HMMscan and HMMsearch are eliminated by parameterization to ensure identical results, just faster.
* TEBinSorter utilizes intelligent parallel workload balancing to optimize the efficiency of searches and ensure near-perfect CPU utilization
* TEBinSorter implements an alternative, optional algorithm for more rapidly searching protein databases (most of what TESorter uses) by pre-screening TEs to identify probable best hits before expensive searches are performed.

Results are stored in a SQLite database, enabling post-hoc filtering and analysis without re-running searches. Run once, filter results at will.

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

## Notes ##

* All results reported by TEBinSorter ultimately emerge from identical HMM alignments to identical sequences with identical locations and hit probability metrics (e.g. E-values) using the same --nobias logic as the original; this guarantees that a hit found by TEBinSorter is bitwise identical to a hit found by TESorter. In the default mode, all hits are always 100% identical.
* In the facets mode, TEBinSorter is finding a relatively small subset of all TESorter hits that are extremely likely to be the single best hit per input sequence. This means that the facets search finds less and stops looking at a sequence very quickly. However, after TESorter filters results, the filtering produces a subset of hits essentially identical (99.98%) to those found by TEBinSorter - TEBinSorter finds these hits directly instead of by post-processing from a completely exhaustive search and skips the cost of the exhaustive search on most sequences.

## Benchmarks

Rice TE library (2,431 sequences), 4 processors.

**Important caveat:** TEsorter baseline times were measured on WSL2 (Windows Subsystem for Linux), which has known I/O and process-spawning overhead compared to native Linux. TEsorter's architecture (many subprocess calls to hmmscan, heavy temp file I/O) is disproportionately affected by this. The true speedup on native Linux may be smaller. TEBinSorter times are also WSL2 but its architecture (in-process pyhmmer, minimal disk I/O) is less sensitive to this penalty. I'm going to measure a true comparsion a bit later - the numbers will not be nearly as impressive as they appear here, unless you too are running the program on WSL.

### Default mode

| Database | Models | TEsorter | TEBinSorter | Speedup |
|----------|--------|----------|-------------|---------|
| TIR | 17 | 26m 41s | 1.0s | ~1,601x |
| LINE | 28 | 56m 45s | 1.8s | ~1,892x |
| SINE | 87 (excl. SINE_SO) | ~30m | 10.4s | ~173x |
| REXdb | 266 | >2h | 18.9s | >381x |
| GyDB | 314 | >2h | 18.7s | >385x |

All 5 databases searched in **50 seconds** total search time, **54 seconds** wall clock (4 processors). SINE_SO excluded by default.

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

Tested on the rice6.9.5.liban TE library (2,431 sequences). Default mode results filtered using TEsorter's exact thresholds (`coverage >= 20`, `evalue <= 1e-3`, `probability >= 0.5`, `normalized score >= 0.1`) and compared at the (sequence, model) pair level. Run with --compat-tesorter-rounding for demonstration purposes.

| Database | TEsorter pairs | TEBinSorter pairs | Common | TEsorter only | Notes |
|----------|---------------|-------------------|--------|---------------|-------|
| TIR | 430 | 430 | 430 | 0 | Perfect match |
| LINE | 806 | 806 | 806 | 0 | Perfect match |
| REXdb | 2,787 | 2,787 | 2,787 | 0 | Perfect match |
| GyDB | 7,496 | 7,496 | 7,496 | 0 | Perfect match |
| SINE | 0 (TEsorter bug) | 731 | -- | -- | TEsorter forces `seq_type='prot'` for a DNA database |

**Full pipeline (HMM + classifier + BLAST pass-2):** 572/584 exact classification match (97.9%) on rice REXdb. 10 diffs from principled score-ordered deconfliction (TEsorter uses arbitrary file-order iteration, which can cause a weaker hit to hold a domain slot against a stronger competitor). 2 misses: borderline sequences whose BLAST pass-2 classification depends on slight differences in the classified target pool.

**Alignment coordinates:** Identical to TEsorter's GFF3 output

**E-values:** Identical (Z set to match hmmscan convention).

**TEsorter rounding bug:** TEsorter rounds `domain_score / model_length` to 2 decimal places before threshold comparison. This results in cases where 9.9 rounds -> 10 and therefore passes a TESorter filter it never should have. TEBinSorter uses full precision by default. `--compat-tesorter-rounding` replicates the old rounding behavior.

### Large-scale benchmarks

133,611 TE consensus sequences (concatenated filtered library), 10 processors, WSL2. The benchmark dataset is available on [Figshare](TODO: add DOI link after upload).

#### Search timing: default vs facet mode

Full pipeline (HMM search + classification + BLAST pass-2) on 4 amino acid databases:

| Database | Models | Default | Facet | Speedup |
|----------|--------|---------|-------|---------|
| REXdb | 266 | 1001s | 488s | **2.1x** |
| GyDB | 314 | 1272s | 926s | **1.4x** |
| LINE | 28 | 66s | 76s | 0.9x |
| TIR | 17 | 55s | 45s | 1.2x |
| **Total wall clock** | | **42min 41s** | **26min 53s** | **1.6x** |

Facet mode per-stage breakdown:

| Database | Facet screen | Verify | Cross-family | Legacy fallback | Total |
|----------|-------------|--------|-------------|----------------|-------|
| REXdb (266) | 277s | 60s | 36s | 97s | 488s |
| GyDB (314) | 546s | 103s | 72s | 185s | 926s |
| LINE (28) | 32s | 8s | 0s | 34s | 76s |
| TIR (17) | 17s | 9s | 0s | 18s | 45s |

LINE and TIR have single-family databases, so no cross-family search is needed.

#### Classification agreement: default vs facet mode

46,464 sequences classified by default mode, 46,324 by facet mode (full pipeline including BLAST pass-2).

| Metric | Count | % of common |
|--------|-------|------------|
| Sequences in both | 45,974 | — |
| Exact match (order+superfamily+clade) | 42,396 | 92.2% |
| Order+superfamily match | 45,143 | 98.2% |
| Order match | 45,688 | 99.4% |
| Order differs | 286 | 0.6% |
| Default only | 490 | — |
| Facet only | 350 | — |

Per-database breakdown:

| Database | Default | Facet | Common | Exact match | Default only | Facet only |
|----------|---------|-------|--------|-------------|-------------|------------|
| REXdb | 38,951 | 38,653 | 38,629 | 33,928 (87.8%) | 322 | 24 |
| GyDB | 37,738 | 37,243 | 37,205 | 32,118 (86.3%) | 533 | 38 |
| LINE | 151 | 147 | 147 | 138 (93.9%) | 4 | 0 |
| TIR | 19,512 | 19,489 | 19,489 | 19,489 (100%) | 23 | 0 |

Most classification differences are at the clade level within the same order and superfamily. These arise from the facet search landing on a different model within the same TE family, which can produce a different clade assignment during deconfliction. TIR classifications are identical between modes.

### SINE_SO: excluded by default

The AnnoSINE database contains SINE_SO (M=4,176), an outlier model 6x larger than the next largest SINE model. Analysis of its contribution:

| Dataset | SINE_SO raw hits | Survive TEsorter filters | % of AnnoSINE compute |
|---------|-----------------|------------------------|----------------------|
| Rice (2,431 seqs) | 2,674 | 0 | 71% |
| 133k sequences | 207,554 | 1 | 71% |

SINE_SO accounts for 71% of AnnoSINE's total M² compute cost but produces effectively zero usable hits after standard filtering. TEBinSorter excludes it by default:

| Dataset | With SINE_SO | Without | Speedup |
|---------|-------------|---------|---------|
| Rice (4 proc) | 39.4s | 10.4s | 3.8x |
| 133k seqs (10 proc) | 38.5min | 10.9min | 3.5x |

Use `--include-sine-so` or `-d sine-so` to search it explicitly.

## Installation

### Dependencies

- Python >= 3.10
- [pyhmmer](https://pyhmmer.readthedocs.io/) >= 0.10
- [pyfastx](https://github.com/lmdu/pyfastx) >= 2.0
- numpy >= 2.0

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
| `{prefix}.db` | SQLite database with all results (hits, classifications, BLAST) |
| `{prefix}.aa` | Six-frame translated amino acid sequences (indexed) |
| `{prefix}.{db}.cls.tsv` | Per-database TE classifications (order, superfamily, clade, completeness) |
| `{prefix}.cls.tsv` | Combined classifications across all databases + BLAST pass-2 |
| `{prefix}.{db}.classifications.tsv` | Facet classifications with confidence tiers (facet mode) |
| `blast_pass2/` | BLAST database and query chunks (temporary) |

## Architecture

### Core modules

- **`pipeline.py`** — Main CLI and search orchestration
- **`search.py`** — HMM search engine with balanced parallelism
- **`sequence.py`** — FASTA ingestion (pyfastx) and six-frame translation (pyhmmer)
- **`hmm.py`** — HMM loading, alphabet detection, optimized profile construction
- **`results.py`** — SQLite persistence with pre-parsed columns (base_seq, strand, frame, domain_type)
- **`deconflict.py`** — Numpy-based hit deconfliction and parameterized filtering

### Facet classification

- **`decompose_hmm.py`** — Standalone sub-HMM decomposition and splicing. Tiered window sizes, configurable overlap. Works with DNA and amino acid HMMs.
- **`model_graph.py`** — Cross-HMM-model similarity graph. Pre-built for all databases. Contains a record of how similar each HMM record is to all the other records in the database.
- **`facet_classify.py`** — Facet screen → verification → cross-family completion → legacy fallback
- **`cross_family.py`** — Targeted search for missing domain families in classified frames

### Classification and post-processing

- **`classifier.py`** — Config-driven TE classification from domain hits. Per-database rules for domain remapping, overlap-aware deconfliction, order/superfamily/clade assignment.
- **`blast_pass2.py`** — BLAST-based pass-2 classification for HMM-unclassified sequences. Parallel chunked BLAST with cross-database target pooling.

# Extended methods

## HMM facets

The concept of an HMM facet is one I've been exploring for other projects but found a use for in this one. By utlizing pyhmmer, we can observe the emission probabilities of an HMM model in a friendly manner in Python. I use these emission probabilities to find conserved subregions of the HMM model which are highly influential in the decisionmaking process of HMMer to actually search a sequence. These high-influence regions are extracted to form a "facet," whose emission probabilities are cloned from the original's over the corresponding window. There are multiple advantages to this approach:

* Short models tend towards specificity. A full-length model will make an effort to compile information about required search effort across an entire sequence, while short models confirm or reject a local region quickly. Facets encoding the same number of amino acids as their parent tend nonetheless to search the same sequence more quickly by parts than the original model.
* Much of the information that the full-length model would glean about the appropriateness of a search from a long view of a sequence is extremely highly correlated with what a good facet will report. A good facet hit almost ensures a good full-length model hit.
* Facets can be intelligently sized so that they exactly pack SIMD lanes (96, 64, 32 amino acids, get consumed mostly in 16-AA sized bites) in the HMMer internals. This reduces low-level CPU waste compared to less politely divisible sizes of sequence.

Are they mathematically correct or statistically sound? I honestly have no idea. But they work.

For completeness, the facet generation code is intentionally a separate module from the remainder of TEBinSorter. It can be used in other projects.

## Parallel strategy

### Legacy (default) search

PyHMMer exposes hmmsearch behavior with two internal, C-level parallelization schemes. These are "queries", corresponding with a search parallelized over the HMM models where each thread picks up one HMM model and searches it against all input sequences, and a "targets" mode where one sequence is searched against each HMM model in parallel. Queries parallelism is inherently more efficient, unless there are only a few models. 

HMM model runtimes scale approximately with their model length M^2. Longer models can therefore dominate runtime despite being only "one" model. AnnoSINE has exactly such a case - the SINE_SO HMM model has a length of ~4100 bp and accounts for ~71% of the total runtime of the SINE model set. When the "queries" parallel mode is used, what ends up happening is that all of the other HMMs in annosine finish rather quickly, and SINE_SO runs single-threaded for usually >10x longer than all the rest of the models took put together. Using the original HMMer, you can provide more threads to a single model - it alleviates the problem somewhat, but HMMer's own internal parallelism is not very efficient.

TEBinSorter wrangles this problem by caclulating the expected runtime cost of each HMM model in a database in advance (cost=M^2), and grouping them into one bin of "close enough in size to run in 'queries' mode" and one bin of "large models that benefit from 'targets' mode. Small models are those whose M^2 is <= the 75th percentile + (2 x IQR) among model costs for that database. Large models are any others.

This all effects low-overhead, near-perfect parallelism in the most efficient available modes. Sequences are reused from the same in-memory object for both searches, so there is essentially no cost to this process.

### Facets search

The facets search is staged:

- Extracted facets from all models are searched with --nobias settings exactly like the full sequences are in the legacy search
- The top facet hits for each sequence are ordered by hit quality, noting the parent HMM model that the facet hit corresponds to.
- Each sequence is searched against its top scoring facet's parent at full length using --nobias. This is exactly the same search that the legacy module does, just targeted instead of exhaustive
- If a sequence lands a full-length model hit, it is "verified" and no further searches are done. If not, the next best scoring facet hit for a different parent HMM model is used. The process repeats until the sequence verifies, or falls out. Most sequences verify against their first facet hit. Almost all verify within 3.
- Because the verification stops at a single best hit, it's possible a sequence may have domains marked by a different family of TEs which are not reached by way of the facet search + verification. The cross-domain search module takes each verified facet hit and searches the same sequence against the parent model associated with the top facet hit of each other TE family. This ensures that the sequence gets both primary and secondary labelling sensitivity both within and across TE families.
- Finally, the remaining sequences which are unclassified by this point, including those with facet hits but which failed any level of verification, are searched in legacy mode.

The expensive parts of HMMer occur exactly when a match is found. What the facets search does is relatively inexpensively route each sequence to exactly and only the HMM models likely to produce its best hits. This saves expensive, and for TEs extremely repetitive, search costs. Most HMM models within the same TE family, e.g. SINEs, will all hit the same sequence at differing strengths. This is exactly when HMMer will burn a lot of runtime, but all of the weaker hits would be discarded later, anyway. It's better (faster) to not do them at all.

The remnant sequences shunted to the legacy search are mostly genuine rejects with no TE signal according to the TESorter databases. These are comparatively inexpensive searches - for example, a rexDB legacy search spends ~92% of its runtime finding the hits within just 71k of about 880k protein frames in our largest testing dataset. They are, mostly, quickly rejected by every model. In a facets search, only the ~6% of sequences missed by the facets have real costs to searching in this way. This is the other major location where a facets search saves time.
