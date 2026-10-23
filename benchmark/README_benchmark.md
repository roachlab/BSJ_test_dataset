# BSJ pipeline benchmark: CIRCexplorer2 + CIRI2 (and your own pipeline)

Runs CIRCexplorer2 and CIRI2 on the eight read libraries in `../dataset/`, scores each output
against the ground truth, and produces a comparison table. The evaluator (`evaluate.py`) is
tool-agnostic, so you can score your own pipeline the same way.

The eight libraries are the full **{PE, SE} x {150, 100} bp x {low, high error}** factorial, so
read length and library layout can be compared independently — hold two axes fixed to isolate
the third (e.g. `pe150_lowerr` vs `se150_lowerr` for paired-ness, `pe150_lowerr` vs
`pe100_lowerr` for read length). Four more libraries (`pe150_lowerr_rev`, `pe150_lowerr_unstr`,
`se150_lowerr_rev`, `se150_lowerr_unstr`) carry the same fragments as their `*150_lowerr` base
but reverse-stranded (dUTP-style, R1 antisense — typical total RNA-seq) or unstranded; the base
libraries are forward-stranded. Trim the `libraries` block in `config.yaml` to run a subset.

## Layout

| File | Purpose |
|------|---------|
| `evaluate.py` | Parse a tool's calls, match to truth, compute metrics. Core of the benchmark. |
| `run_circexplorer2.sh` | STAR align (chimeric) -> CIRCexplorer2 parse + annotate -> evaluate, per library. |
| `run_ciri2.sh` | BWA-MEM -> CIRI2 -> evaluate, per library. |
| `run_benchmark.sh` | Run both tools, then summarise. |
| `combine_results.py` | Merge all `eval_*.json` into one table + `results/summary.tsv`. |
| `run_bsj_finder.sh` | Drive the real **bsj_finder** pipeline (detection only) per library -> convert -> evaluate. |
| `bsj_finder_to_bed.py` | Convert bsj_finder's per-read `bsj_cor.tsv.gz` into an aggregated BED6 for scoring. |

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
python3 evaluate.py --tool bed --calls mypipe_pe150_lowerr.bed \
    --dataset ../dataset --library pe150_lowerr --out-prefix results/mypipeline/eval_pe150_lowerr
python3 combine_results.py results
```

## Benchmarking bsj_finder

bsj_finder is the pipeline this dataset was built to validate, so it is wired in as a
first-class tool. It is itself a Snakemake pipeline living one level up from this submodule
(the `BSJ_test_dataset` repo is a submodule of `bsj_finder`). The benchmark drives it to its
**detection target only** (`bsj_results/<lib>.bsj_cor.tsv.gz`), skipping the mmseqs2 dedup and
bcftools polish steps, which refine circle *sequences* and do not change breakpoint
*coordinates* — so they cannot affect recall/precision/specificity.

bsj_finder emits **one row per read** that spans a back-splice, with genomic breakpoints
(`s_start`/`s_end`) already ordered low->high and on the same 0-based convention as
`truth_bsj.tsv`. `bsj_finder_to_bed.py` clusters reads whose breakpoints agree within
`--merge-tol` bp (default 5) into one junction, sums the supporting reads (→ the
quantification metric), and writes a BED6 that `evaluate.py --tool bed` scores like any other
tool.

```bash
# needs snakemake + bsj_finder's own deps (--use-conda builds the per-rule envs)
cd benchmark
BSJ_FINDER=.. THREADS=8 ./run_bsj_finder.sh ../dataset results/bsj_finder
python3 combine_results.py results
```

Or via Snakemake (bsj_finder is in the `tools:` list in `config.yaml`):

```bash
snakemake --cores 8 results/bsj_finder/eval_pe150_lowerr.json   # one library
snakemake --cores 8                                             # everything, all tools
```

Settings in `config.yaml`: `bsj_finder_dir` (repo root, relative to `benchmark/`; default `..`)
and `bsj_finder_merge_tol` (breakpoint clustering tolerance). Coverage note: the two
**J2** isoforms share identical breakpoints and are scored as one truth junction (`BSJ_J2`),
so 48/48 is reachable by any coordinate-based caller, bsj_finder included.

**Presets as two tools.** `bsj_finder_variants` in `config.yaml` maps tool names to extra
bsj_finder settings: `bsj_finder` (its default `sensitive` preset, `-k 11 -w 5`) and
`bsj_finder_fast` (`preset=fast`, `-k 21 -w 11`). Each is benchmarked as its own tool under
`results/<name>/`, with its own minimap2 index, and both are in `tools:` by default. The presets
need bsj_finder `cc4da4f` or later: an older checkout silently ignores `preset=`, so the index
step stops with an error if bsj_finder does not report the preset it applied. Point
`bsj_finder_dir` at a checkout with presets; the default, the `bsj_finder` submodule, is pinned
to `cc4da4f`, which has them.

**Failed runs are recorded, not fatal.** If bsj_finder exits non-zero for a library (e.g.
strandedness auto-detection refusing on a small, noisy library), the benchmark writes an empty
call set and `results/<tool>/<lib>.status` = `failed: <reason>`; `evaluate.py --status` carries
it into the JSON and `summary.tsv` gains a `status` column (`ok` otherwise). bsj_finder's own
log for each library is `results/<tool>/work/<lib>/bsj_finder.log`.

## What gets matched

A call matches a true BSJ when chrom is equal and **both** breakpoints fall within a tolerance
window of the truth coordinates. The window is `max(--base-tol, breakpoint_tolerance_bp)`, where
`breakpoint_tolerance_bp` comes from `truth_bsj.tsv` (non-zero for the fuzzy-breakpoint H and
indel-junction L cases). `--base-tol` (default 2) absorbs 0/1-based and end-convention
differences between tools. Each truth BSJ and each call is matched at most once.

## Metrics (per tool x library)

Positives = the 48 true BSJ loci; negatives = the 39 decoy loci.

- **Recall** = TP / 48 — fraction of true BSJs detected. Also broken out **per category**
  (A-M) and **per coverage bin** (1-2, 3-9, 10-49, 50+ junction reads) so you can read the
  sensitivity floor directly.
- **Precision** = TP / (TP + FP). FP is split into `FP_on_decoy` (a call landing on a decoy
  locus — the meaningful false positives) and `FP_offtarget` (calls landing nowhere near truth
  or a decoy).
- **F1** = harmonic mean of precision and recall — the single headline number.
- **Specificity** = TN / 39 over the decoy set (decoys correctly not called).
- **Accuracy** = (TP + TN) / (48 + 39) over the positive+decoy locus universe. Reported for
  completeness; F1 + specificity are more informative for detection because off-target FPs are
  not part of this denominator.
- **FDR** = FP / (TP + FP).
- **Quantification** = Pearson + Spearman correlation of detected vs `expected_junction_reads`
  across matched calls — how well junction-read counts track the truth.
- **Strand** — matching is strand-blind; the strand each tool declares for a matched call is
  then compared to `truth_bsj.tsv`: `strand_correct` / `strand_wrong` / `strand_unknown` (`.`),
  `strand_accuracy` = correct / TP, `recall_strand` = correct-strand TP / 48, and
  `strand_accuracy_by_truth_strand` (`+` vs `-`; the truth set is 42 `+` / 6 `-`, so always
  answering `+` scores 0.875). `strandedness` (fwd / rev / unstranded) is taken from
  `../dataset/libraries.tsv`. Read the `*_rev` / `*_unstr` runs against their base library:
  a caller that reports read-alignment strand is right on fwd, wrong on every call on rev.
- **missed** / **decoys_called** / **strand_wrong_ids** lists (in the JSON) name exactly which
  loci failed, for drill-down.

## Runtime and memory

Each tool rule wraps its main command in GNU `time` and writes
`results/timing/<tool>/<lib>.<rule>.tsv` (wall s, user s, system s, max RSS KB). One-off index
builds are `shared.<rule>.tsv`. `plot_results.py` turns these into `plots/runtime.png` (wall time,
CPU time and peak memory per tool) and `results/runtime_summary.tsv`, **for the scaled libraries
only** (`scale` > 1 in `libraries.tsv`): on the small libraries every step takes seconds and
mostly measures start-up.

- **Use the scaled libraries.** On the 3,568-read libraries every step takes seconds and is
  dominated by start-up, so they are left out. `python3 generate_dataset.py --scaled` writes
  `{pe,se}150_lowerr_x{100,1000,10000}` (same truth; ~357k / 3.6M / 36M reads or pairs, ~1.4 GB
  in total) to `dataset/scaled/`, rendered as independently seeded 100x chunks on up to 32
  processes (about 7 min). `python3 verify.py --scaled` checks them. They are git-ignored.
  Only the x100 pair is listed in `config.yaml` (as `optional`); add the others to run them.
- **Genome size.** `python3 generate_dataset.py --padded-genomes` writes
  `dataset/genomes/genome_1gbp.fa`: the three real contigs unchanged plus 100 random 10 Mbp
  filler contigs (~45% GC), so reads, truth and annotation stay valid. The workflow still uses
  `python3 verify.py --padded-genomes` checks it. Git-ignored.
- **Running against another genome / a subset of libraries.** `genome=`, `results=` and `only=`
  select the reference, a separate results tree and the libraries, e.g.
  `snakemake --cores 8 --config threads=8 genome=../dataset/genomes/genome_1gbp.fa
  results=results_genome_1gbp 'only=["pe150_lowerr_x1000","se150_lowerr_x1000"]'`. The x1000 /
  x10000 libraries are `scaled_only` and run only when named in `only=`. bsj_finder's minimap2
  index is built once per tree (`bsj_finder_index`) and linked into each library run, so like the
  STAR and BWA indexes it is timed separately. `bsj_finder_config="minimap2_k=21 minimap2_w=11"`
  passes extra settings to bsj_finder's own config (index and detection), e.g. its `sr`-preset
  seeding for large genomes; use a separate `results=` tree to keep both runs.
- **Scaling plots.** `python3 plot_scaling.py --tree results:1.5e6 --tree
  results_genome_1gbp:1e9 --runtime-genome results_genome_1gbp --out results/plots` writes
  `runtime_by_library_size` (wall and CPU time vs reads, log-log, on the 1 Gbp genome) and
  `memory_by_genome_size` (peak memory vs genome size for the x100 runs and the index builds),
  plus `scaling_summary.tsv`. bsj_finder is shown for both of its seeding presets: the base trees
  use `sensitive` (its default), and `--variant` adds `fast` as a dashed line from one results tree
  per genome size, with `--rename` labelling the base line:

  ```bash
  python3 plot_scaling.py --tree results:1.5e6 --tree results_genome_1gbp:1e9 \
      --runtime-genome results_genome_1gbp --dataset ../dataset --out results/plots \
      --rename "bsj_finder=bsj_finder (sensitive)" \
      --variant "bsj_finder (fast)|bsj_finder|results_k21:1.5e6,results_genome_1gbp_k21:1e9"
  ```

  The `fast` trees come from runs with `results=results_k21` / `results_genome_1gbp_k21` and
  `bsj_finder_config="preset=fast"` (bsj_finder with presets) or, equivalently,
  `bsj_finder_config="minimap2_k=21 minimap2_w=11"`.
- **Run one job at a time with the same threads for every tool**, e.g.
  `snakemake --cores 8 --config threads=8`, or the timings measure contention.
- **Container:** the current `bsj-bench.sif` lacks GNU `time`; bind the host's:
  `apptainer exec -B /usr/bin/time:/usr/bin/time ...`. The recipes now install it, so a rebuilt
  image does not need this. Override the path with `--config time_bin=...`.
- CPU and peak RSS cover child processes; RSS is the **largest single process**, not a sum.
- bsj_finder's fastp step used to run `--dedup --dup_calc_accuracy 4`, which allocates a fixed
  ~8 GB hash regardless of input and, without UMIs, removed independent junction reads that
  share a sequence (recovery fell to 0.22-0.41 of expected reads on the scaled libraries). With
  dedup dropped, its peak is ~1.2-1.4 GB on the scaled libraries; most of the remaining ~1 GB is
  fastp's duplication-rate report (`--dont_eval_duplication` would remove it).

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
`results/circexplorer2/`, `results/ciri2/`, `results/bsj_finder/`, `results/summary.tsv`, `results/plots/*.png`.

Edit `config.yaml` to change threads, the breakpoint tolerance, the CIRI2 path, the read
libraries, or which `tools` to run. To benchmark only one tool, drop the other from `tools:`
and request the targets you want (e.g. `snakemake --cores 8 results/ciri2/eval_pe150_lowerr.json`).
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

### `std::bad_alloc` during `singularity build`

This is an out-of-memory / out-of-disk error in the build temp area, not a recipe bug.
Singularity's temp and cache default to `/tmp`, which is often small or RAM-backed (tmpfs);
packing the squashfs or running the conda solve then exhausts it. Fixes:

```bash
# (a) give the build a roomy, disk-backed temp + cache, then rebuild
export APPTAINER_TMPDIR=/data/$USER/tmp  APPTAINER_CACHEDIR=/data/$USER/cache
#     older Singularity: SINGULARITY_TMPDIR / SINGULARITY_CACHEDIR
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
singularity build bsj-bench.sif bsj-bench.def

# (b) low on RAM? build with Docker (heavy work there) and just convert to SIF
sudo docker build -t bsj-bench .
sudo singularity build bsj-bench.sif docker-daemon://bsj-bench:latest
```
