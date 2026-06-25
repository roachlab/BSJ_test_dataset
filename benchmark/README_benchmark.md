# BSJ pipeline benchmark: CIRCexplorer2 + CIRI2 (and your own pipeline)

Runs CIRCexplorer2 and CIRI2 on the four read libraries in `../dataset/`, scores each output
against the ground truth, and produces a comparison table. The evaluator (`evaluate.py`) is
tool-agnostic, so you can score your own pipeline the same way.

## Layout

| File | Purpose |
|------|---------|
| `evaluate.py` | Parse a tool's calls, match to truth, compute metrics. Core of the benchmark. |
| `run_circexplorer2.sh` | STAR align (chimeric) -> CIRCexplorer2 parse + annotate -> evaluate, per library. |
| `run_ciri2.sh` | BWA-MEM -> CIRI2 -> evaluate, per library. |
| `run_benchmark.sh` | Run both tools, then summarise. |
| `combine_results.py` | Merge all `eval_*.json` into one table + `results/summary.tsv`. |

## Environment

The aligners/callers are not bundled. Easiest is the pinned conda env
(`environment.yml`: Python 3.11, Snakemake 8, STAR, BWA, samtools, CIRCexplorer2,
gtfToGenePred, perl):

```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate bsj-bench
# CIRI2 is distributed as a Perl script (not on bioconda). Download CIRI2.pl and:
export CIRI2=/path/to/CIRI2.pl
```

Or skip local setup entirely and use Docker (see below) — it bakes in every tool,
**including CIRI2**.

## Run

```bash
cd benchmark
./run_benchmark.sh ../dataset          # runs both tools + prints the summary
# or individually:
THREADS=8 ./run_circexplorer2.sh ../dataset results/circexplorer2
THREADS=8 ./run_ciri2.sh        ../dataset results/ciri2
python3 combine_results.py results
```

## Scoring your own pipeline

Emit a BED6 per library — `chrom  start(0-based)  end  name  junction_reads  strand` — where
`start`/`end` are the circle's left/right genomic breakpoints, then:

```bash
python3 evaluate.py --tool bed --calls mypipe_pe_lowerr.bed \
    --dataset ../dataset --library pe_lowerr --out-prefix results/mypipeline/eval_pe_lowerr
python3 combine_results.py results
```

## What gets matched

A call matches a true BSJ when chrom is equal and **both** breakpoints fall within a tolerance
window of the truth coordinates. The window is `max(--base-tol, breakpoint_tolerance_bp)`, where
`breakpoint_tolerance_bp` comes from `truth_bsj.tsv` (non-zero for the fuzzy-breakpoint H and
indel-junction L cases). `--base-tol` (default 2) absorbs 0/1-based and end-convention
differences between tools. Each truth BSJ and each call is matched at most once.

## Metrics (per tool x library)

Positives = the 49 true BSJ loci; negatives = the 39 decoy loci.

- **Recall** = TP / 49 — fraction of true BSJs detected. Also broken out **per category**
  (A-M) and **per coverage bin** (1-2, 3-9, 10-49, 50+ junction reads) so you can read the
  sensitivity floor directly.
- **Precision** = TP / (TP + FP). FP is split into `FP_on_decoy` (a call landing on a decoy
  locus — the meaningful false positives) and `FP_offtarget` (calls landing nowhere near truth
  or a decoy).
- **F1** = harmonic mean of precision and recall — the single headline number.
- **Specificity** = TN / 39 over the decoy set (decoys correctly not called).
- **Accuracy** = (TP + TN) / (49 + 39) over the positive+decoy locus universe. Reported for
  completeness; F1 + specificity are more informative for detection because off-target FPs are
  not part of this denominator.
- **FDR** = FP / (TP + FP).
- **Quantification** = Pearson + Spearman correlation of detected vs `expected_junction_reads`
  across matched calls — how well junction-read counts track the truth.
- **missed** / **decoys_called** lists (in the JSON) name exactly which loci failed, for
  drill-down.

## Expected behaviour / gotchas

- CIRCexplorer2 `annotate` reports only **annotation-supported** circles, so it will miss the
  un-annotated (C), intronic (D), EIciRNA (E), and intergenic/antisense (M) cases by design.
  The script therefore evaluates **both** the de-novo `back_spliced_junction.bed`
  (`*_denovo`) and the annotated `circularRNA_known.txt` (`*_known`); compare the two.
- CIRI2 uses the GTF only to label circle type; it calls de-novo, so it can recover the
  un-annotated cases.
- The non-canonical (G), fuzzy (H), and indel (L) cases stress each tool's motif and
  breakpoint handling — expect tool-specific differences there.
- Very-low-coverage circles (1-2 reads) sit at the single-read detection limit; differences in
  the 1-2 coverage bin are where pipelines most diverge.

## Snakemake workflow (recommended)

`Snakefile` + `config.yaml` run everything as a DAG: references -> alignment -> calling ->
evaluation -> `summary.tsv` + plots, with per-library parallelism and automatic resume.

```bash
conda activate bsj-bench          # provides Snakemake 8 + all tools (see environment.yml)
export CIRI2=/path/to/CIRI2.pl    # if CIRI2.pl is not on PATH

cd benchmark
snakemake -n --cores 8            # dry run: print the DAG
snakemake --cores 8               # run the whole benchmark
```

Outputs land in `results/` exactly as with the shell scripts:
`results/circexplorer2/`, `results/ciri2/`, `results/summary.tsv`, `results/plots/*.png`.

Edit `config.yaml` to change threads, the breakpoint tolerance, the CIRI2 path, the read
libraries, or which `tools` to run. To benchmark only one tool, drop the other from `tools:`
and request the targets you want (e.g. `snakemake --cores 8 results/ciri2/eval_pe_lowerr.json`).
The shell scripts (`run_*.sh`) remain as a Snakemake-free alternative.

## Singularity / Apptainer (recommended on HPC / no-root machines)

Singularity runs as **your** user (no root daemon) and auto-mounts your working directory, so
the volume-permission problems that bite Docker on shared machines simply don't occur — outputs
are written to `benchmark/results/` owned by you. The definition file is `bsj-bench.def` at the
repository root; it bakes in every tool, **including CIRI2**.

```bash
# build (from the repository root, so %files can see benchmark/environment.yml)
sudo singularity build bsj-bench.sif bsj-bench.def
#   unprivileged alternative if fakeroot is configured:
#   singularity build --fakeroot bsj-bench.sif bsj-bench.def
#   or convert an existing Docker image:
#   singularity build bsj-bench.sif docker-daemon://bsj-bench:latest

# run from the repository root (CWD is auto-mounted and writable as you)
singularity run  bsj-bench.sif                # snakemake --cores 4 in benchmark/
singularity run  bsj-bench.sif --cores 8      # pass extra snakemake args
singularity exec bsj-bench.sif bash -c "cd benchmark && snakemake --cores 8"
```

Use `apptainer` in place of `singularity` if that is your command. Building still needs root or
`--fakeroot`; everything after the build is unprivileged.

## Docker (alternative; zero-setup, includes CIRI2)

`Dockerfile` (at the repository root) builds an image with every tool plus Snakemake 8 and
fetches CIRI2 at build time.

```bash
# from the repository root (build context must be the repo root)
docker build -t bsj-bench .

# run the whole benchmark; results land in ./benchmark/results on the host
docker run --rm -v "$PWD/benchmark/results:/work/benchmark/results" bsj-bench

# more cores / custom command
docker run --rm -v "$PWD/benchmark/results:/work/benchmark/results" \
    bsj-bench bash -c "cd benchmark && snakemake --cores 8"
```

If SourceForge is unreachable at build time, pass another CIRI2 source with
`--build-arg CIRI2_URL=<url>`, or mount `CIRI2.pl` into the container and set `CIRI2`.

### Volume permissions

The container runs as the non-root `mambauser` (uid 1000). If you bind-mount a host
`benchmark/results` that is root-owned (e.g. Docker created it under `sudo`), the in-container
user cannot write to it and you get `permission denied` creating `results/star_index`. Pick one:

```bash
# (a) run as root (simplest when you already use sudo); outputs come back root-owned
sudo docker run --rm --user root \
    -v "$PWD/benchmark/results:/work/benchmark/results" bsj-bench
sudo chown -R $USER benchmark/results          # reclaim ownership afterwards

# (b) pre-create the dir for the container uid; outputs owned by you, no --user needed
mkdir -p benchmark/results && sudo chown -R 1000:1000 benchmark/results
sudo docker run --rm -v "$PWD/benchmark/results:/work/benchmark/results" bsj-bench
```

Or skip the mount entirely and copy results out afterwards with `docker cp`.


## Troubleshooting

### STAR "illegal instruction" on older CPUs (no AVX2)

The bioconda STAR binary is compiled with AVX2 and crashes with SIGILL on CPUs that lack it
(e.g. Intel Ivy Bridge / Xeon E5 v2). bwa, samtools and the callers are unaffected. The
Singularity (`bsj-bench.def`) and Docker (`Dockerfile`) recipes now rebuild STAR from source
with `-msse4.2`, which runs on those CPUs (negligible speed cost on this small genome) — just
rebuild the image:

```bash
singularity build bsj-bench.sif bsj-bench.def        # or: docker build -t bsj-bench .
# even older CPU? drop to plain SSE:
#   docker build --build-arg STAR_SIMD=sse -t bsj-bench .
#   (Singularity: change STAR_SIMD in bsj-bench.def)
```

For a **bare conda env** (no container), replace the AVX2 STAR in place:

```bash
conda activate bsj-bench
conda install -y -c conda-forge make gxx_linux-64 zlib    # build tools, if missing
VER=2.7.11b
curl -L "https://github.com/alexdobin/STAR/archive/refs/tags/${VER}.tar.gz" | tar xz
make -C "STAR-${VER}/source" STAR CXXFLAGS_SIMD="-msse4.2"
cp "STAR-${VER}/source/STAR" "$CONDA_PREFIX/bin/STAR"     # overwrite the AVX2 build
STAR --version                                            # sanity check
```
