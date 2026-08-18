# TEsorter2
 
Fast, divergence-robust classification of transposable elements.
 
TEsorter2 is a reimplementation of [TEsorter](https://github.com/zhangrengang/TEsorter) that
keeps its classification semantics while introducing three major improvements:
 
1. **Speed**: HMMER's limited parallelism is replaced by [pyHMMER](https://pyhmmer.readthedocs.io/)
   with workload-aware load balancing, plus an optional facet pre-screen and a parallelized
   second-pass BLAST.
2. **Sensitivity on degraded copies**: an optional [BATH](https://github.com/TravisWheelerLab/BATH)
   engine performs frameshift-aware translated search directly against nucleotide sequence.
3. **Reproducibility**: clade assignment uses a score-weighted vote instead of a count-based vote
   whose ties were broken by internal data-structure ordering.
It also adds multi-database reconciliation in a single run and a BATH-based genome mode.

## Installation

### conda (recommended)

Installs the external binaries (HMMER, BLAST+) and TEsorter2 in one step:

```bash
git clone https://github.com/KGerhardt/TEsorter2.git
cd TEsorter2
conda env create -f environment.yml
conda activate tesorter2
```

### pip

Requires Python >= 3.9. HMMER and BLAST+ must already be on `PATH`:

```bash
pip install git+https://github.com/KGerhardt/TEsorter2.git
```

### Databases

The HMM databases (REXdb, GyDB2, LINE, TIR, AnnoSINE) ship inside the package, as they do in
TEsorter, so there is no download step and no configuration: `tesorter2 sequences input.fasta` works
straight after install.

> **GyDB2 is shipped modified.** Its 2009-era profiles carried no `COMPO` record, which HMMER and
> pyHMMER tolerate but stricter HMMER3 readers (notably `nail`) reject outright. The bundled copy
> adds `COMPO` to all 314 profiles, derived from each model's own match emissions by HMMER's
> `p7_hmm_SetComposition()` — no probability was altered and no model was rebuilt. Verified to give
> byte-identical hits and classifications. Full notice, rationale and reproduction steps in
> [`tesorter2/database/LICENSES.md`](tesorter2/database/LICENSES.md).

To use a custom collection of HMM databases instead, point TEsorter2 at its directory:

```bash
tesorter2 sequences input.fasta --db-dir /path/to/db   # or: export TESORTER2_DB=/path/to/db
```

Individual databases can also be passed by path: `-d /path/to/custom.hmm`.

### BATH (required for `--genome`, optional for `--bath`)

[BATH](https://github.com/TravisWheelerLab/BATH) is not packaged on conda, and despite what its
`INSTALL` file suggests, it publishes no release artifacts — the repo carries tags but zero
releases. A source build is the only route. Easel is a separate clone, not a submodule:

```bash
git clone https://github.com/TravisWheelerLab/BATH.git
cd BATH
git clone -b BATH https://github.com/TravisWheelerLab/easel
autoconf && ./configure && make
```

Binaries land in `BATH/src`. Put that on `PATH`, or set `BATH_BIN_DIR` to it. Verified against
BATH 2.0 (`v2.0.0-rc4`); needs only gcc/make/autoconf.

Genome mode searches protein profiles with BATH exclusively, so BATH is a hard requirement whenever
`--genome` is used with an amino-acid database. A DNA-only genome run (`-d sine`) uses nhmmer and
does not need BATH.

### nail (optional, for cascade stages)

[nail](https://github.com/TravisWheelerLab/nail) seeds with MMseqs2 and approximates HMMER's
Forward/Backward. It is only needed when a `--stages` list includes `nail`:

```bash
cargo install nail        # needs MMseqs2 on PATH; nail shells out to it for seeding
```

It reads HMMER3 `.hmm` files directly, so there is no conversion step. Amino-acid only: it cannot
search DNA profile databases and has no role in genome mode. Verified against nail 0.7.1 with
MMseqs2 18.8cc5c.

## Two modes

The mode is a verb, because what the input *is* decides nearly everything downstream — which
engines can run, whether the input is windowed, whether there is a per-element classification at
all:

```bash
tesorter2 sequences lib.fa      # pre-extracted TE sequences -> per-element classifications
tesorter2 genome    asm.fa      # whole assembly -> domain-level GFF3 + summary
```

The older flag form (`tesorter2 lib.fa`, `tesorter2 asm.fa --genome`) still works and warns.

## Choosing an engine
 
| | Element mode (pre-extracted TEs) | Genome mode (assembly) |
|---|---|---|
| **Intact / low-divergence sequence** | **pyHMMER** (default) : fastest | **BATH** (only engine) |
| **Degraded / frameshifted copies** | **BATH** (`--bath`) : slower than pyHMMER, more sensitive | **BATH** (only engine) : both faster *and* more sensitive |

Genome mode offers no engine choice: BATH is both faster and more sensitive on long sequences, so
the pyHMMER genome path was removed rather than kept as a slower, less sensitive alternative.

---
 
## Search modes
 
### Default mode (pyHMMER)
 
Single-pass `--nobias` search against all models via pyHMMER, in-process: `hmmsearch` for
amino-acid databases, `nhmmer` for DNA databases (`sine`, `sine-so`).
 
- Model-cost-aware parallel load balancing chooses between pyHMMER's `queries` and `targets`
  parallelization per model bin, reaching near-full CPU utilization (see
  [Parallel strategy](#parallel-strategy)).
- Results are written to a SQLite database, so filtering and re-analysis do not require re-running
  the search.
  
### DNA databases (nhmmer)

DNA profile databases are searched with `nhmmer`, not `hmmsearch`. `hmmsearch` only scores the
strand it is handed, so it misses every minus-strand copy, and pyHMMER rejects sequences over
100k residues outright. `nhmmer` scans both strands in one pass and handles long targets.

On 5 Mb of rice against `sine`, nhmmer records hits on both strands (15,992 `+` / 15,821 `-`)
where hmmsearch records no strand at all, and classifies 50 windows against hmmsearch's 37. The
35 windows both engines classify are almost all `+` or unstranded; every one of the 15 windows
only nhmmer recovers is on the minus strand, so the gain is the strand hmmsearch cannot see
rather than a looser threshold.

`--dna-engine hmmsearch` restores the old single-strand behaviour for comparison.

**Known limitation — hit filters on the DNA path are inherited from the protein path.** DNA hits
are filtered with the same thresholds as protein domains (coverage ≥ 20%, E-value ≤ 1e-3,
accuracy ≥ 0.5, normalized score ≥ 0.1). Two consequences on the rice/`sine` fixture:

- **Two of the four filters never discriminate.** Accuracy and coverage reject nothing: the same
  50 windows are classified whether the cutoffs are at their defaults or at zero. Only E-value
  and normalized score bind.
- **The E-value cutoff is not comparable between engines.** nhmmer scores against a long-target
  search space spanning both strands, hmmsearch against a per-sequence protein-style one, so the
  same alignment gets very different E-values — for `SHANSINE_MT` on one rice window, 5.5e-05
  under hmmsearch versus 0.0014 under nhmmer. A single `1e-3` cutoff is therefore stricter for
  nhmmer than for hmmsearch, and the two windows hmmsearch classifies that nhmmer does not are
  both boundary misses of this kind, not detection failures — nhmmer finds hundreds of raw hits
  in each.

The practical effect is bounded: dropping the E-value filter altogether raises the count from 50
to 59, so no threshold choice recovers much more. SINE-specific filters are not implemented; use
`--dna-engine hmmsearch` if you need the old behaviour for comparison.

### Facet pre-screen (`--stages facet`)

Pre-screens amino-acid databases with spliced sub-HMMs ("facets") to route each sequence only to
the models likely to produce its best hit:

1. **Facet screen**: tiered sub-HMMs (96 → 64 → 48 → 32 aa) searched against all six translated
   frames.
2. **Targeted verification**: top facet hit per domain family verified with a full-model
   `--nobias` search.
3. **Cross-family completion**: verified frames searched for missing domain families.
4. **Legacy fallback**: frames with no facet signal get a full search.

Faster than a plain `hmmer` stage on AA databases, with 99.8% post-filter recall. Available as a
cascade engine only — `--stages facet` or `--stages facet,bath`. It is mutually exclusive with
`hmmer`: facet *is* hmmsearch with a pre-screen in front, and step 4 already runs the full search
on whatever the screen missed.

DNA databases always use nhmmer; DNA facets do not repay their overhead.

### BATH mode (`--bath`)
 
Replaces pyHMMER with `bathsearch --fs` for amino-acid databases. BATH aligns protein pHMMs
directly against raw nucleotide input, allowing the alignment to change frame at indels and read
through stop codons under a frameshift penalty, so a domain interrupted by a frameshift is
recovered as a single hit with its true extent. No six-frame translation is performed.
 
- HMMER3 databases are converted on first use (`bathconvert`) and cached as `{db}.bath.hmm`.
- Hits are normalized into the same internal schema as the HMMER path, so classification,
  reconciliation and BLAST pass-2 are unchanged.
- BATH tblout reports no per-residue posterior probability, so `acc` is set to `1.0` and the
  minimum-accuracy filter is a no-op for BATH hits.
- Minus-strand coordinates are normalized to ascending order, with strand encoded in the target
  suffix, matching the HMMER convention.
- Implied by `--genome`, where BATH is the only protein engine.

### Genome mode (`--genome`)
 
Treats the input as whole-genome sequence rather than pre-extracted elements: detects TE protein
domains throughout, classifies **each domain individually**, resolves overlapping features, and
emits a domain-level GFF3 plus a summary table. It does not produce a per-element `.cls.tsv` and
does not run BLAST pass-2, matching TEsorter's `-genome` behaviour.
 
Protein domains are found with **BATH only**. `bathsearch --fs` runs directly on the nucleotide
windows; the tblout already reports nucleotide coordinates and strand, so no six-frame translation
or coordinate back-mapping happens at all. BATH additionally streams long targets in ~0.25 Mb
blocks overlapped by the maximum expected hit length, reconciling boundary duplicates internally.
`--bath` is therefore implied and has no effect here.

DNA databases (`sine`) are searched with nhmmer against the **raw, un-windowed** input, so
non-coding elements are annotated independently of the protein path's window sizing. A DNA-only run
is supported and needs no BATH; windowing is skipped entirely and no `cut.fa` is written.

Overlaps are resolved within each feature type, so a BATH protein-domain feature (`CDS`) and an
nhmmer element feature (e.g. `SINE_element`) at the same locus coexist as independent layers rather
than evicting one another.

Window size: `--win-size` (default `1e6`), `--win-ovl` (default `1e5`).
Requires at least one database.
 
### Multiple databases
 
**Every bundled database is searched by default.** Omitting `-d` searches rexdb, gydb, line, tir
and sine together; pass `-d rexdb,gydb` to restrict the set. Each database classifies
independently and reconciliation resolves them afterwards, so more databases means more evidence
rather than more ambiguity — at proportionally more compute.

The three reconciliation stages:
 
1. **Independent classification**: each database classifies every element on its own, emitting a
   per-database `{prefix}.{db}.cls.tsv` in native TEsorter format.
2. **Name harmonization**: per-database calls are projected onto a unified taxonomy that collapses
   superfamily synonyms (e.g. `Pao` → `Bel-Pao`) and unifies clade names. Lineages with no
   established equivalent are kept distinct to avoid spurious agreement.
3. **Scope-aware hierarchical reconciliation**: a hierarchical vote at Order → Superfamily →
   Clade. At each level candidate labels are weighted by the summed normalized domain score per
   database, and only entries consistent with the winning label advance. The superfamily level is
   *scope-aware*: a database may only elect a superfamily it models with at least 2 clades.
Every per-database call is retained in the `SecondaryHits` column as
`db:order/superfamily/clade=score`, in descending order of evidence — always on the database's
**native** names, so the audit trail survives harmonization. Use `--compat-tesorter-output`
for the original 7-column format.

> **Compatibility contract.** The per-database `{prefix}.{db}.cls.tsv` is **never harmonized** and
> always carries exactly TEsorter's original 7 columns. Downstream pipelines — EDTA reads
> `*.{db}db.cls.tsv` positionally as `(id, Order, Superfamily)` — depend on this, so harmonized
> names and any added columns belong in the combined `{prefix}.cls.tsv` only. Note that
> harmonization rewrites `Order` and `Superfamily` too, not just `Clade` (`Pao`→`Bel-Pao`,
> `pararetrovirus`→`LTR`/`Caulimoviridae`), which is precisely why it must not reach the per-database
> file. If you ever feed the *combined* file to EDTA, extend the `%lib` hash in its
> `cleanup_misclas.pl` first — it knows `pararetrovirus` but not `Caulimoviridae` or `Bel-Pao`.

Harmonization is gated on at least two databases voting on an element: a single-database run keeps
native names throughout, so `-d rexdb` alone is unaffected. The two tables driving stages 2 and 3
are `tesorter2/database/clade_harmonization.tsv` (mapping) and `clade_scope.tsv` (which superfamily
each database can resolve at lineage level). The REXdb↔GyDB lineage equivalences — `Ale`/`Retrofit`,
`Ivana`/`Oryco`, `Tekay`/`Del`, `SIRE`, `Tork`, `Reina`, `CRM`, `Galadriel`, `Athila` — are asserted
from [Neumann et al. 2019](https://doi.org/10.1186/s13100-018-0144-1) (*Mob DNA* 10:1), not fitted
to any dataset. Both tables are read from the resolved `--db-dir` first and fall back to the
packaged copies, so a custom database collection still harmonizes.

Where one database resolves *below* the level both can express, the combined file adds a `Lineage`
column rather than forking the clade name. REXdb splits the Tat group into `Ogre`/`Retand`/`TatI-III`
while GyDB has a single `tat`: both are reported as `Clade=Tat`, and REXdb's finer call lands in
`Lineage`. Lineage is chosen among the databases that actually resolve one, so the finer call
survives even when the coarser database wins the score vote. (`Tatius` is not part of the group —
the REXdb HMM places it at `OTA/Tatius`, a sibling of `Tat`.)
 
### Sequence Ontology

`Order/Superfamily/Clade` is TEsorter's vocabulary, not a standard one. Every classification is
also resolved to a [Sequence Ontology](http://www.sequenceontology.org) term, so results are
comparable with other annotation tools:

- `{prefix}.cls.tsv` gains `SO_name` and `SO_ID` columns (e.g. `Copia_LTR_retrotransposon`,
  `SO:0002264`).
- Genome-mode GFF3 features carry `Ontology_term=SO:...` plus `so_name=`. The feature type stays
  `CDS`: these features are protein domains, not elements, so typing one as
  `Gypsy_LTR_retrotransposon` would assert the domain *is* the retrotransposon.

The mapping authority is [EDTA's `TE_Sequence_Ontology.txt`](https://github.com/oushujun/EDTA/blob/master/bin/TE_Sequence_Ontology.txt),
bundled in `tesorter2/data/`. All 59 Order/Superfamily labels the bundled databases can emit
resolve to a specific SO term; nothing falls back to the generic `repeat_region`.

For lineages with no SO term of their own, EDTA files a descriptive name under a generic
accession (`CR1_LINE_retrotransposon` is not a real SO term; its `SO:0000194` is
`LINE_element`'s). `SO_name` keeps EDTA's name for interoperability, while anything written as an
ontology term resolves to the real one.

`--compat-tesorter-output` suppresses the SO columns, keeping the original 7-column format.

### Clade voting
 
Within each database the winning clade is chosen by a **score-weighted vote**: each domain
contributes a length-normalized score (`dom_score / model_len`) to its clade, and the highest
summed score wins:
 
$$\hat{c} = \arg\max_{c} \sum_{d \in D_c} \frac{\text{score}_d}{\text{len}_d}$$
 
Normalizing by model length controls for the tendency of longer profiles to accumulate higher raw
scores, making domains comparable. Ties fall through to the existing mixture/completeness rules;
Order and Superfamily are inherited from the winning clade.
 
This replaces TEsorter's raw domain-count plurality, which breaks ties by position in an internal
collection rather than by any biological signal, causing sibling-clade swaps (Reina↔Tekay,
Ale↔Alesia) to flip depending on search engine. Pass `--compat-tesorter-voting` to restore the
original behaviour.
 
---

## CLI reference
 
```
tesorter2 <sequence> [options]
```
 
| Flag | Default | Description |
|---|---|---|
| `sequence` | — | Input FASTA (TE library, or genome with `--genome`) |
| `-d`, `--database` | all bundled | Comma-separated database aliases or paths |
| `--max-search` | off | Force the full database set even alongside an explicit `-d` |
| `-o`, `--outdir` | `{input}.TEsorter2` | Output directory |
| `--db-dir` | bundled | Directory holding the HMM databases (see Installation) |
| `--dna-engine` | `nhmmer` | Engine for DNA databases (`nhmmer` or `hmmsearch`) |
| `--prefix` | input basename | Output file prefix |
| `-p`, `--processors` | `4` | Processors |
| `--bath` | off | Frameshift-aware BATH engine (AA databases only; implied by `--genome`) |
| `--stages` | off | Staged multi-engine cascade over AA databases, e.g. `nail,hmmer,bath` |
| `--mask-stops` | off | Translate stop codons as `X` instead of `*`, letting profiles align through premature stops |
| `--genome` | off | Genome mode: domain-level annotation + GFF3 (BATH required) |
| `--win-size` | `1e6` | Genome mode window size |
| `--win-ovl` | `1e5` | Genome mode window overlap |
| `--emit-bath` | off | Emit routed FASTA partitions for BATH to `{outdir}/BATHwater/` |
| `--include-sine-so` | off | Include the SINE_SO model in AnnoSINE searches |
| `--compat-tesorter-voting` | off | Raw domain-count clade vote (TEsorter behaviour) |
| `--compat-tesorter-rounding` | off | Round normalized scores to 2 dp before filtering (replicates a TEsorter bug) |
| `--compat-tesorter-output` | off | Emit combined `.cls.tsv` in TEsorter's 7-column format |

---
 
## Output files
 
| File | Description |
|---|---|
| `{prefix}.db` | SQLite database with all hits, classifications and BLAST results |
| `{prefix}.aa` | Six-frame translated amino-acid sequences (indexed; HMMER path only) |
| `{prefix}.{db}.cls.tsv` | Per-database classifications (order, superfamily, clade, completeness) |
| `{prefix}.cls.tsv` | Combined classifications across databases + BLAST pass-2 (+ `SecondaryHits`, `SO_name`, `SO_ID`) |
| `{prefix}.dom.gff3` | Genome mode: classified TE protein-domain features |
| `{prefix}.dom.fna` | Genome mode: nucleotide sequences of each classified feature |
| `{prefix}.genome.summary.tsv` | Genome mode: Order/Superfamily/Clade tallies |
| `blast_pass2/` | BLAST database and query chunks (temporary) |
| `cut.fa` | Genome mode: windowed genome, written only when a protein database is searched (temporary) |
 
---
 
## Benchmarks
 
All runs on an ANVIL CPU node (Rosen Center for Advanced Computing, Purdue University;
AMD EPYC 7763, 256 GB RAM), 16 threads, mean ± SD over 3 replicates. BATH in frameshift-aware
mode (`--fs`).
 
### Element mode
 
108,318 RepBase elements (LTR, TIR, LINE), each searched against its order-specific database
(GyDB, REXdb-pnas, REXdb-line).
 
| Pipeline | Engine | Time | Speedup vs TEsorter |
|---|---|---|---|
| TEsorter | HMMER (`hmmscan`) | 2,932 s | 1.0× |
| **TEsorter2** | **pyHMMER** | **543 s** | **5.4×** |
| TEsorter2 | BATH | 983 s | 3.0× |
 
### Genome mode
 
Complete *Oryza sativa* genome (~375 Mb) against REXdb.
 
| Pipeline | Engine | Time | Speedup vs TEsorter |
|---|---|---|---|
| TEsorter | HMMER (`hmmscan`) | 2,127 s | 1.0× |
| TEsorter2 | pyHMMER *(path since removed)* | 988 s | 2.2× |
| **TEsorter2** | **BATH** | **431 s** | **4.9×** |
 
BATH is 2.3× faster than pyHMMER here because it avoids the six-frame translation and coordinate
back-mapping that dominate the HMMER path on long sequences. Being both faster and more sensitive
on long targets is why the pyHMMER genome path was removed; the row is retained as the measurement
that motivated the decision.

---
 
## Output files
 
| File | Description |
|---|---|
| `{prefix}.db` | SQLite database with all hits, classifications and BLAST results |
| `{prefix}.aa` | Six-frame translated amino-acid sequences (indexed; HMMER path only) |
| `{prefix}.{db}.cls.tsv` | Per-database classifications (order, superfamily, clade, completeness) |
| `{prefix}.cls.tsv` | Combined classifications across databases + BLAST pass-2 (+ `SecondaryHits`, `SO_name`, `SO_ID`) |
| `{prefix}.dom.gff3` | Genome mode: classified TE protein-domain features |
| `{prefix}.dom.fna` | Genome mode: nucleotide sequences of each classified feature |
| `{prefix}.genome.summary.tsv` | Genome mode: Order/Superfamily/Clade tallies |
| `blast_pass2/` | BLAST database and query chunks (temporary) |
| `cut.fa` | Genome mode: windowed genome, written only when a protein database is searched (temporary) |
 
---
 
TEsorter2 | BATH | 983 s | 3.0× |

---

## TEsorter compatibility
 
`tesorter2-compat` provides a drop-in CLI with TEsorter's original argument names and
defaults (including count-based clade voting), for substituting TEsorter inside existing pipelines:
 
```bash
tesorter2-compat input.fasta -db rexdb -p 16 -pre out
```
 
Supported: `-db/--hmm-database`, `--db-hmm`, `-st/--seq-type`, `-pre/--prefix`, `-p/--processors`,
`-tmp/--tmp-dir`, `-cov/--min-coverage`, `-eval/--max-evalue`, `-prob/--min-probability`,
`-score/--min-score`, `-dp2/--disable-pass2`, `-nolib/--no-library`, `-norc/--no-reverse`,
`-nocln/--no-cleanup`.
 
---
 
## Architecture
 
### Core
 
- **`pipeline.py`** — CLI and search orchestration
- **`search.py`** — HMM search engine with balanced parallelism
- **`sequence.py`** — FASTA ingestion (pyfastx) and six-frame translation (pyhmmer)
- **`hmm.py`** — HMM loading, alphabet detection, optimized profile construction
- **`bath_search.py`** — BATH engine: conversion, `bathsearch --fs` invocation, hit normalization
- **`genome.py`** — Genome mode: windowing, per-domain classification, overlap resolution, GFF3/summary
- **`results.py`** — SQLite persistence with pre-parsed columns (base_seq, strand, frame, domain_type)
- **`deconflict.py`** — numpy-based hit deconfliction and parameterized filtering
### Facet classification
 
- **`decompose_hmm.py`** — standalone sub-HMM decomposition and splicing (tiered windows,
  configurable overlap; works for DNA and AA HMMs)
- **`model_graph.py`** — cross-model similarity graph, precomputed for bundled databases
- **`facet_classify.py`** — screen → verify → cross-family completion → legacy fallback
- **`cross_family.py`** — targeted search for missing domain families in classified frames
### Classification and post-processing
 
- **`classifier.py`** — config-driven classification from domain hits; per-database domain
  remapping, overlap-aware deconfliction, order/superfamily/clade assignment
- **`blast_pass2.py`** — parallel chunked BLAST pass-2 with cross-database target pooling
- **`tesorter_output.py`**, **`emit.py`**, **`id_registry.py`** — output formatting and identifier bookkeeping
---
 
## Extended methods
 
### HMM facets
 
pyHMMER exposes an HMM's emission probabilities in Python. TEsorter2 uses them to locate conserved
subregions that most influence HMMER's decision-making, and extracts each such region into a
"facet": a complete, self-contained HMM whose emission and transition probabilities are cloned from
its parent over the corresponding window.
 
Facets are sized at 96, 64, 48 or 32 amino acids, using the longest size the parent supports;
models of ≤32 aa are used as-is. Windows are chosen to maximize the summed emission probability
over the parent and may overlap by at most 33% of facet length.
 
For a parent HMM of length *M* and a query of length *N*, one `hmmsearch` is a dynamic-programming
matrix of size *M × N*. A facet search replaces *M* with a facet size *F* (*F* < *M*), yielding an
*F × N* slice with a proportional reduction in work. Three properties make this profitable:
 
- **Short models are more decisive.** A full-length model accumulates information about search
  effort across the whole sequence; a short model confirms or rejects a local region quickly.
- **Facet scores predict full-length scores.** A good facet hit almost always implies a good
  full-length hit, making facets a fast approximation of the final score.
- **Facet sizes pack SIMD lanes.** 96/64/32 aa are consumed in 16-aa bites by HMMER's internals,
  reducing low-level CPU waste relative to less divisible sizes.
Facets do not reproduce their parent's domain detections exactly, so a full-length verification is
still required for correctness. The acceleration comes from *skipping* verifications: the TE HMMs
within each database are highly redundant (either wholly, as in single-type databases, or in
subcollections, as in REXdb and GyDB), so most models in a cluster hit the same sequence at
differing strengths and all but the best are discarded downstream anyway — yet detecting a weak hit
costs exactly as much as detecting a strong one. TEsorter2 ships precomputed all-vs-all similarity
graphs (obtained by searching each model's consensus against every other model in the database)
that identify these clusters. Facet hits order the parent searches within each cluster from highest
to lowest, and verification stops as soon as a high-scoring parent is confirmed — typically 1–2
parent searches per cluster per sequence instead of dozens.
 
In effect, the deconfliction that would otherwise happen *after* an exhaustive search is moved
*before* it, at the cost of a cheap approximate screen. Post-filter agreement with the exhaustive
search is 99.98% at hit level and 99.8% at family recall (4 misses of 2,090 on rice REXdb).
 
The facet generation code is intentionally a separate module and is reusable in other projects.
 
### Parallel strategy
 
**Default search.** pyHMMER exposes two C-level parallelization schemes: `queries`, where each
thread takes one HMM and searches it against all sequences, and `targets`, where one sequence is
searched against models in parallel. `queries` is inherently more efficient unless there are few
models.
 
HMM runtime scales roughly with *M²*, so a single long model can dominate. AnnoSINE is the
pathological case: `SINE_SO` (M≈4,100) accounts for ~71% of the model set's runtime, and under
`queries` parallelism every other model finishes quickly while `SINE_SO` runs single-threaded for
>10× longer than all the rest combined.
 
TEsorter2 precomputes each model's expected cost (*M²*) and bins them: **small** models
(cost ≤ 75th percentile + 2×IQR for that database) run in `queries` mode; **large** models run in
`targets` mode. Sequences are reused from the same in-memory object across both searches, so the
split is essentially free, yielding near-full CPU utilization in the most efficient mode available
for each model class.
 
**Facet search.** Staged: facets from all models are searched with `--nobias` exactly as in the
default search; top facet hits per sequence are ordered by quality, noting each facet's parent;
each sequence is searched against its best facet's parent at full length. If it verifies, no
further searches run; otherwise the next-best facet hit for a *different* parent is tried. Most
sequences verify on the first attempt; almost all within three. Because verification stops at a
single best hit, a cross-family pass then searches each verified sequence against the top facet
hit's parent for every *other* TE family, preserving both primary and secondary labelling
sensitivity. Sequences still unclassified fall back to the legacy search.
 
The expensive part of HMMER is finding a match. Facets route each sequence only to the models
likely to produce its best hits, skipping the redundant weak hits that would be discarded anyway.
The leftovers sent to the legacy search are mostly genuine rejects with no TE signal, and are
cheap: a REXdb legacy search spends ~92% of its runtime on the ~71k of ~880k protein frames that
actually contain hits.
 
---



## License

GPL-3.0-or-later. See [LICENSE](https://github.com/KGerhardt/TEsorter2/blob/master/LICENSE).

The bundled HMM databases (REXdb, GyDB2, AnnoSINE, Kapitonov LINE, Yuan & Wessler TIR) are
third-party data with their own upstream licenses — CC BY 4.0 (REXdb), Creative Commons
Attribution (GyDB2), MIT (AnnoSINE), and redistribution via TEsorter/GPL-3.0 (Kapitonov LINE,
Yuan & Wessler TIR). TEsorter2's GPL-3.0 does **not** extend to them. Per-database licenses,
sources, and required citations are in
[`tesorter2/database/LICENSES.md`](https://github.com/KGerhardt/TEsorter2/blob/master/tesorter2/database/LICENSES.md); cite the databases you run
against.
