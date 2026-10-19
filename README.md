# BSJ Test Dataset

A small, fully-labelled synthetic benchmark for **back-splice-junction (BSJ / circRNA)
detection pipelines**. Every read's origin is recorded, so any caller can be scored for
**sensitivity** (recall of true BSJs), **specificity** (rejection of engineered decoys),
and **quantification** (junction-read counts). Reproducible from a single seed
(`generate_dataset.py`, `SEED=20260624`).

The dataset was built to validate [`bsj_finder`](https://github.com/roachlab/bsj_finder) —
included here as a git submodule — and is wired up against **CIRCexplorer2** and **CIRI2**
as reference callers in `benchmark/`.

---

## Repository layout

```
.
├── generate_dataset.py      # generates everything in ./dataset/ (seeded, ~700 lines)
├── verify.py                # sanity checks: coords in-bounds, motifs match genome, read counts
├── BSJ_test_dataset_plan.md # original design document
├── Dockerfile               # one-image benchmark env (STAR/BWA/CIRCexplorer2/CIRI2/Snakemake)
├── bsj-bench.def            # Singularity/Apptainer definition of the same env
├── bsj-bench.sif            # (built image; not normally committed)
├── dataset/                 # the dataset (see below)
│   ├── README.md            # dataset file manifest + coordinate conventions + scoring recipe
│   ├── genome.fa            # 1.5 Mbp toy genome: chrT1 (700 kb), chrT2 (600 kb), chrP (200 kb)
│   ├── annotation.gtf       # gene models — deliberately omits some features (see below)
│   ├── reads_*.fastq.gz     # 12 read libraries: 8-way factorial + 4 strandedness variants
│   ├── libraries.tsv        # per-library layout / length / error / strandedness
│   ├── truth_bsj.tsv        # 48 true BSJs (ground truth)
│   ├── truth_bsj.bed        # BED6 of the two breakpoints per circle
│   ├── truth_bsj_sequences.fa  # junction windows (±50 bp) + full circle sequences
│   ├── truth_negatives.tsv  # 39 decoy loci — expected result: no circle call
│   └── read_origin.tsv      # per-read provenance for scoring
├── benchmark/               # tool-agnostic benchmark harness
│   ├── README_benchmark.md  # detailed benchmark docs (env setup, gotchas, troubleshooting)
│   ├── evaluate.py          # core scorer: parse calls → match truth → metrics (tool-agnostic)
│   ├── combine_results.py   # merge eval_*.json → results/summary.tsv
│   ├── plot_results.py      # plots: recall vs coverage, per-category, P/R/F1, quantification
│   ├── run_circexplorer2.sh # STAR (chimeric) → CIRCexplorer2 parse+annotate → evaluate
│   ├── run_ciri2.sh         # BWA-MEM → CIRI2 → evaluate
│   ├── run_bsj_finder.sh    # drive the real bsj_finder pipeline (detection only) → evaluate
│   ├── bsj_finder_to_bed.py # convert bsj_finder per-read bsj_cor.tsv.gz → aggregated BED6
│   ├── Snakefile            # Snakemake DAG: index → align → call → evaluate → summary + plots
│   ├── config.yaml          # threads, tolerances, CIRI2 path, tool list, library definitions
│   └── environment.yml      # pinned conda env (bsj-bench)
└── bsj_finder/              # git submodule — the pipeline under test (its own Snakemake pipeline)
```

---

## The dataset

### Reference genome

Three contigs, **1,500,000 bp** total; ~45% GC random background with deliberately placed
exons, introns, repeats, and junction motifs:

| Contig | Length | Role |
|--------|--------|------|
| `chrT1` | 700 kb | Most positive-case genes |
| `chrT2` | 600 kb | Remaining positives + most decoys |
| `chrP` | 200 kb | Paralog / multi-mapping decoys (near-identical copies vs `chrT2`) |

### Ground truth: 48 true BSJs (`truth_bsj.tsv`)

Each positive category is instantiated at **multiple independent loci**, with junction
coverage spread across orders of magnitude (down to 1–2 reads) so the sensitivity floor can
be mapped per case type:

| Cat | Category | n | What it tests |
|-----|----------|---|---------------|
| A | Canonical multi-exon | 5 | Positive control / typical circle (GT/AG) |
| B | Single-exon circle | 4 | Self-circularization of one exon |
| C | Un-annotated exon | 4 | Circle uses an exon **absent from the GTF** |
| D | Intronic ciRNA | 4 | Circle contained in an intron, no exon overlap |
| E | EIciRNA (exon-intron) | 3 | Retained intron in the circle body |
| F | Minus-strand | 4 | Strand handling / motif on − strand |
| G | Non-canonical motif | 3 | Motif filters (AT/AC etc. instead of GT/AG) |
| H | Fuzzy breakpoint | 4 | Microhomology at the breakpoint; truth carries a tolerance window |
| I | Sub-read-length (rolling) | 4 | Circle shorter than the read; reads cross the BSJ ≥2× |
| J | Alternative back-splicing | 3 | Multiple circles sharing/alternating breakpoints |
| K | Nested / overlapping circles | 3 | Two circles from one gene nesting inside each other |
| L | Indel at junction | 3 | Inserted/deleted bases exactly at the back-splice |
| M | Intergenic (2) + antisense (2) | 4 | No annotation / strand assignment |

Motif breakdown: 37 canonical GT/AG, 11 non-canonical (incl. all of G, most of D/L).
Strand: 42 `+`, 6 `−` (F and the two M antisense circles).

**Coverage spread:** junction-read targets range from 200 (A1) down to 1 (D4), including
several 1–2 read cases at the single-read detection limit.

**Matching tolerances:** `breakpoint_tolerance_bp` is non-zero only for H (3–6 bp,
microhomology) and L (2–4 bp, indel) — those loci accept calls anywhere in the window.

> **48 vs 49:** the original plan listed 49 because category J2 had two body isoforms
> (`BSJ_J2_1`/`BSJ_J2_2`) sharing identical breakpoints. They are collapsed into a single
> junction `BSJ_J2`, because this benchmark scores *back-splice-junction discovery*, not
> isoform resolution. A coordinate-based matcher can therefore match at most one of the two
> J2 isoforms.

### Ground truth: 39 decoy loci (`truth_negatives.tsv`)

Engineered look-alikes that a correct caller must **not** call. Critically, the false
breakpoints sit on **canonical GT/AG motifs**, so the trap cannot be escaped by motif
filtering alone:

| Decoy | n | The trap |
|-------|---|----------|
| DEC1 | 4 | Tandem duplication; rotated reads look back-spliced |
| DEC2 | 5 | Normally spliced linear gene; forward SJs must not be miscalled as BSJ |
| DEC3 | 4 | Near-identical paralogs (`chrP`/`chrT2`) → multi-mapping chimeras |
| DEC4 | 3 | Inverted repeat / palindrome mimicking a chimeric junction |
| DEC5 | 3 | Readthrough / trans-splice: chimeric-looking but linear |
| DEC6 | 3 | Homopolymer / low-complexity soft-clip artifacts |
| DEC7 | 4 | *Joint test*: spliced gene with a tandem repeat inside it — keep real forward SJs, reject the repeat-induced pseudo-BSJ |
| DEC8 | 4 | *Joint test*: gene with near-duplicate exons — keep real splicing, reject the within-gene pseudo-BSJ |
| DEC9 | 4 | RT template-switching chimera (library artifact, not genomic) |
| DEC10 | 3 | SNP near a real breakpoint; must not create a spurious isoform/shift |
| DEC11 | 2 | Very high-expression linear gene swamping the background (dynamic-range stress) |

DEC7/DEC8 are **joint** tests: the same locus emits both `host_linear` forward-SJ reads
(expected: retained) and `negative` pseudo-BSJ reads (expected: no circle call);
`read_origin.tsv` records both roles.

### Read libraries

Full factorial **{PE, SE} × {150, 100} bp × {low, high error} = 8 libraries**. Layout and
read length vary *independently*, so a PE-vs-SE difference can be attributed to the
mate-pair signal rather than read length (and vice versa). Every library contains the same
**3,568 reads from the same underlying truth** — only layout, length, and error differ
(PE: 330 bp fragment, R1/R2 files).

| Library | Layout | Length | Error model |
|---------|--------|--------|-------------|
| `pe150_lowerr` / `pe150_higherr` | paired | 2×150 | ~0.2% sub / ~2% sub + ~0.3% indel |
| `pe100_lowerr` / `pe100_higherr` | paired | 2×100 | same |
| `se150_lowerr` / `se150_higherr` | single | 150 | same |
| `se100_lowerr` / `se100_higherr` | single | 100 | same |

**Deliberate read-length case:** `BSJ_I3` is a 120 bp circle — between the two read
lengths — so it *rolls* (reads wrap the loop, crossing the BSJ multiple times) in the 150 bp
libraries but behaves as an ordinary circle in the 100 bp ones. `BSJ_I1/I2/I4` (90/70/60 bp)
roll in every library. Comparing I3 across the two lengths measures read-length dependence
directly.

To isolate one axis, hold the others fixed — e.g. `pe150_lowerr` vs `se150_lowerr` isolates
paired-ness; `pe150_lowerr` vs `pe100_lowerr` isolates read length.

#### Strandedness

All eight libraries above are **forward-stranded**: R1 (or the single read) is the
transcript-sense sequence (ligation-type protocols; featureCounts `-s 1`, HISAT2 `FR`,
Salmon `ISF`/`SF`). Most total-RNA-seq (dUTP, e.g. TruSeq Stranded Total RNA) is
**reverse-stranded**: R1 is antisense (featureCounts `-s 2`, HISAT2 `RF`, Salmon `ISR`/`SR`).
Four extra libraries re-emit the **exact same fragments** as `pe150_lowerr` / `se150_lowerr`
with only the strand protocol changed, so any difference against the base library is due to
strandedness alone:

| Library | Derived from | Strandedness | How |
|---------|--------------|--------------|-----|
| `pe150_lowerr_rev` | `pe150_lowerr` | reverse | R1/R2 swapped |
| `pe150_lowerr_unstr` | `pe150_lowerr` | unstranded | each fragment fwd or swapped (50/50) |
| `se150_lowerr_rev` | `se150_lowerr` | reverse | read reverse-complemented |
| `se150_lowerr_unstr` | `se150_lowerr` | unstranded | each read as-is or reverse-complemented (50/50) |

`libraries.tsv` lists every library's layout, read length, error model and strandedness.
The evaluator scores the **strand each caller declares** for every detected BSJ against
`truth_bsj.tsv`. That catches a caller that finds the junction but reports the read-alignment
strand instead of the transcript strand, which is invisible in a forward-stranded library and
flips every call in a reverse-stranded one.

### Annotation omissions (by design)

`annotation.gtf` deliberately **omits** category-C circle exons, all intronic (D) /
EIciRNA (E) features, and the intergenic/antisense (M) circles — so annotation-dependent
tools (e.g. CIRCexplorer2 `annotate`) miss them, while de-novo callers can recover them.
This is intentional and lets you separate annotation-guided from de-novo performance.

### Per-read provenance (`read_origin.tsv`)

`read_id, library, origin, read_class, spans_junction, expected_role, note`

- `read_class`: `circle_junction` (10,688), `circle_body` (3,360), `chimera` (3,608),
  `linear` (10,888)
- `expected_role`: `positive` (true circle), `host_linear` (real linear transcript),
  `negative` (decoy — should not yield a circle call)

### Coordinate conventions

- All positions in `truth_bsj.tsv` are **0-based**.
- `acceptor_pos` = genomic coordinate of the circle's 5′ start; `donor_pos` = the circle's
  3′ end. The back-splice joins **donor → acceptor**.
- `+` strand: `acceptor < donor`; `−` strand: `acceptor > donor`.
- `truth_bsj.bed` reports `[min(bp), max(bp)+1)` (0-based, BED standard).
- `splice_motif` is derived directly from `genome.fa` at the two breakpoints (strand-aware),
  so it always matches the sequence.

---

## Regenerating & verifying

```bash
python3 generate_dataset.py   # writes everything into ./dataset/ (seeded, deterministic)
python3 verify.py             # asserts coords in-bounds, motifs match genome, read counts
```

---

## Benchmarking

The evaluator (`benchmark/evaluate.py`) is **tool-agnostic**: it normalises a tool's output
to `(chrom, left_bp, right_bp, junction_reads)`, matches calls against `truth_bsj.tsv` with
a per-junction tolerance window, and splits FPs into on-decoy vs off-target.

### Matching rule

A call matches a true BSJ when the chrom is equal and **both** breakpoints fall within
`max(--base-tol, breakpoint_tolerance_bp)` of the truth. `--base-tol` (default 2) absorbs
0/1-based and end-convention differences. Both truth and calls are first
**left-aligned** against the genome to canonicalise microhomology-shifted breakpoints, so a
call reported at the other end of a microhomology run still matches. Each truth BSJ and
each call matches at most once.

### Metrics (per tool × library)

- **Recall** = TP / 48; also per category (A–M) and per coverage bin (1-2 / 3-9 / 10-49 / 50+)
- **Precision** = TP / (TP + FP); FP split into `FP_on_decoy` (call on a decoy locus — the
  meaningful FPs) and `FP_offtarget`
- **F1** — headline single number
- **Specificity** = TN / 39 over the decoy set
- **Accuracy** = (TP + TN) / 87 (positive + decoy universe)
- **FDR** = FP / (TP + FP)
- **Quantification** = Pearson + Spearman correlation of detected vs
  `expected_junction_reads` across matched calls
- **Strand** (scored separately, so recall/precision are unchanged): for matched calls,
  `strand_correct` / `strand_wrong` / `strand_unknown` (no strand reported),
  `strand_accuracy` = correct / TP, `recall_strand` = correct-strand TP / 48, and
  `strand_accuracy_by_truth_strand` (+ vs −). The split by truth strand matters because the
  truth set is 42 `+` / 6 `−`, so a caller that always says `+` still gets 0.875.
- The JSON output also lists exactly which loci were `missed` / `decoys_called` /
  `strand_wrong_ids`

### Environment (pick one)

1. **Conda:** `conda env create -f benchmark/environment.yml` (env `bsj-bench`; STAR, BWA,
   samtools, CIRCexplorer2, Snakemake 8, gtfToGenePred). CIRI2 is **not** on bioconda —
   download `CIRI2.pl` and `export CIRI2=/path/to/CIRI2.pl`.
2. **Docker (zero-setup, includes CIRI2):** `docker build -t bsj-bench .` then
   `docker run --rm -v "$PWD/benchmark/results:/work/benchmark/results" bsj-bench`
   (see `benchmark/README_benchmark.md` for volume-permission fixes).
3. **Singularity/Apptainer:** `singularity build bsj-bench.sif bsj-bench.def` — runs as your
   user (no root daemon), preferred on HPC.

Both container recipes rebuild **STAR from source without AVX2** (SSE4.2), so the image runs
on older CPUs where the AVX2 bioconda build SIGILLs.

### Running

**Snakemake (recommended):**

```bash
cd benchmark
snakemake -n --cores 8        # dry run
snakemake --cores 8           # whole DAG: references → align → call → evaluate → summary + plots
snakemake --cores 8 results/ciri2/eval_pe150_lowerr.json   # or a single target
```

**Shell scripts (Snakemake-free):**

```bash
cd benchmark
./run_benchmark.sh ../dataset                    # CIRCexplorer2 + CIRI2 + summary
THREADS=8 ./run_circexplorer2.sh ../dataset results/circexplorer2
THREADS=8 ./run_ciri2.sh        ../dataset results/ciri2
python3 combine_results.py results
```

**bsj_finder** (the pipeline this dataset was built to validate):

```bash
cd benchmark
BSJ_FINDER=.. THREADS=8 ./run_bsj_finder.sh ../dataset results/bsj_finder
python3 combine_results.py results
```

It drives the real bsj_finder Snakemake pipeline (via its submodule, `--use-conda` builds the
per-rule envs) to the **detection target only** (`bsj_results/<lib>.bsj_cor.tsv.gz`),
skipping the mmseqs2 dedup and bcftools polish steps — those refine circle *sequences* and
cannot change breakpoint coordinates, so they cannot affect recall/precision/specificity.
`bsj_finder_to_bed.py` clusters per-read breakpoints within `--merge-tol` bp (12 in
`config.yaml` / `run_bsj_finder.sh`) into one junction call, sums supporting reads, and
emits a BED6 scored by `evaluate.py --tool bed`. Clustering is strand-blind, and each call
gets the majority `s_strand` of its reads, so the strand is only right if bsj_finder's
per-read strand is the transcript strand (i.e. corrected for library strandedness).

**Your own pipeline:** emit a BED6 per library —
`chrom  start(0-based)  end  name  junction_reads  strand` — where start/end are the
circle's left/right genomic breakpoints and strand is the circle's **transcript** strand
(`.` if unknown), then:

```bash
python3 benchmark/evaluate.py --tool bed --tool-name mypipe \
    --calls mypipe_pe150_lowerr.bed --dataset dataset \
    --library pe150_lowerr --out-prefix benchmark/results/mypipe/eval_pe150_lowerr
python3 benchmark/combine_results.py benchmark/results
python3 benchmark/plot_results.py --results benchmark/results --out benchmark/results/plots
```

### Outputs

`benchmark/results/`: per-tool `eval_<library>.json` (full metrics + missed/decoy lists),
`summary.tsv` (one row per tool × library), and `plots/`
(`recall_vs_coverage.png`, `per_category_recall.png`, `precision_recall_f1.png`,
`quantification.png`, `strand_accuracy.png`).

### Current results (2026-06-25, Docker; 150 bp libraries)

| Tool | Library | Calls | TP | Precision | Recall | F1 | Specificity | Quant (Spearman) |
|------|---------|-------|----|-----------|--------|----|-------------|------------------|
| circexplorer2 (denovo) | pe_lowerr | 34 | 32 | 0.941 | 0.653 | 0.771 | 1.000 | — |
| circexplorer2 (denovo) | pe_higherr | 23 | 23 | 1.000 | 0.469 | 0.639 | 1.000 | — |
| circexplorer2 (denovo) | se_lowerr | 48 | 45 | 0.938 | 0.918 | 0.928 | 1.000 | — |
| circexplorer2 (denovo) | se_higherr | 63 | 45 | 0.714 | 0.918 | 0.804 | 0.949 | — |
| circexplorer2 (known) | pe_lowerr | 28 | 28 | 1.000 | 0.571 | 0.727 | 1.000 | 0.995 |
| circexplorer2 (known) | pe_higherr | 21 | 21 | 1.000 | 0.429 | 0.600 | 1.000 | 0.792 |
| circexplorer2 (known) | se_lowerr | 35 | 35 | 1.000 | 0.714 | 0.833 | 1.000 | 0.994 |
| circexplorer2 (known) | se_higherr | 34 | 34 | 1.000 | 0.694 | 0.819 | 1.000 | 0.993 |
| ciri2 | pe_lowerr | 35 | 35 | 1.000 | 0.714 | 0.833 | 1.000 | 0.996 |
| ciri2 | pe_higherr | 34 | 34 | 1.000 | 0.694 | 0.819 | 1.000 | 0.996 |
| ciri2 | se_lowerr | 35 | 35 | 1.000 | 0.714 | 0.833 | 1.000 | 0.998 |
| ciri2 | se_higherr | 33 | 33 | 1.000 | 0.673 | 0.805 | 1.000 | 0.993 |

Takeaways so far: CIRCexplorer2's **de-novo** mode is layout-sensitive (SE ≫ PE here,
high-error PE collapses to 0.47 recall); its **known** mode is annotation-capped (as
expected, since it reports only annotation-supported circles). CIRI2 is precision 1.0 with
~0.67–0.71 recall and near-perfect quantification. The 100 bp libraries have not been run
yet, and **bsj_finder has not been benchmarked yet** — those are the open items.

### Expected behaviour / gotchas

- CIRCexplorer2 `annotate` misses C/D/E/M by design (annotation omissions); the benchmark
  scores **both** the de-novo BED (`*_denovo`) and the annotated output (`*_known`).
- CIRI2 calls de-novo (GTF only labels circle type), so it can recover the un-annotated cases.
- G (non-canonical), H (fuzzy), and L (indel) stress motif/breakpoint handling —
  tool-specific differences are expected there.
- 1–2 read circles sit at the single-read detection limit; that coverage bin is where
  pipelines diverge most.
- 48/48 recall is the coordinate-matching ceiling (J2 counts once).
- Strand on the forward-stranded 150 bp libraries (from the 2026-06-25 outputs; any length or
  layout): CIRI2 and CIRCexplorer2 `known` take the strand from the host gene, so they report
  the two antisense circles (`BSJ_M_as1/2`) as `+`. CIRCexplorer2's de-novo BED reports `+` for
  all six minus-strand circles (F1–F4, M_as1/2).
  The `*_rev` / `*_unstr` libraries have not been run yet.

---

## Submodule

`bsj_finder/` is a git submodule (`https://github.com/roachlab/bsj_finder.git`). When
cloning this repo, fetch it with `git submodule update --init`; the benchmark's bsj_finder
rules and `run_bsj_finder.sh` expect the bsj_finder repo root to contain a `Snakefile`,
`workflow/`, and `config/`.
