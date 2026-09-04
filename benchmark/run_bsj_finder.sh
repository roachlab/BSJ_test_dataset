#!/usr/bin/env bash
# Run the REAL bsj_finder pipeline (detection only) on the tiny BSJ test dataset and
# evaluate against ground truth. bsj_finder is a Snakemake pipeline, so we drive it to
# its detection target (bsj_results/<lib>.bsj_cor.tsv.gz), convert that per-read table
# to an aggregated BED6, and score it with the tool-agnostic evaluator.
#
# Requires: snakemake (with bsj_finder's own deps, e.g. --use-conda builds the per-rule
# conda envs: fastp/trimnami, seqtk, minimap2, biopython). See README_benchmark.md.
#
#   BSJ_FINDER   path to the bsj_finder repo root (default: ../.. — the parent of this
#                submodule). Must contain Snakefile + workflow/ + config/.
#   THREADS      threads for snakemake (default 4)
#   MERGE_TOL    breakpoint clustering tolerance for the BED converter (default 12)
#   SMK_ARGS     extra args passed to snakemake (default: --use-conda)
set -euo pipefail

DATASET="${1:-../dataset}"
OUT="${2:-results/bsj_finder}"
THREADS="${THREADS:-4}"
MERGE_TOL="${MERGE_TOL:-12}"
SMK_ARGS="${SMK_ARGS---use-conda}"

# Resolve to absolute paths (bsj_finder runs in its own working dir).
DATASET="$(cd "$DATASET" && pwd)"
GENOME="$DATASET/genome.fa"
BENCH="$(pwd)"
mkdir -p "$OUT"
OUT_ABS="$(cd "$OUT" && pwd)"

BSJ_FINDER="${BSJ_FINDER:-../..}"
BSJ_FINDER="$(cd "$BSJ_FINDER" && pwd)"

command -v snakemake >/dev/null 2>&1 || { echo "ERROR: 'snakemake' not found in PATH (see README_benchmark.md)"; exit 1; }
[ -f "$BSJ_FINDER/Snakefile" ] || { echo "ERROR: no Snakefile at BSJ_FINDER='$BSJ_FINDER' (set \$BSJ_FINDER to the bsj_finder repo root)"; exit 1; }
[ -f "$GENOME" ] || { echo "ERROR: genome not found: $GENOME"; exit 1; }
echo "bsj_finder repo : $BSJ_FINDER"
echo "dataset         : $DATASET"

run_lib () {
  local lib="$1"; shift
  local reads=("$@")
  echo ">>> bsj_finder :: $lib"

  # Per-library read dir with deterministic names so metasnek assigns sample == <lib>.
  local rdir="$OUT_ABS/reads/$lib"
  local work="$OUT_ABS/work/$lib"
  rm -rf "$rdir"; mkdir -p "$rdir" "$work"
  if [ "${#reads[@]}" -eq 2 ]; then
    ln -sf "${reads[0]}" "$rdir/${lib}_R1.fastq.gz"
    ln -sf "${reads[1]}" "$rdir/${lib}_R2.fastq.gz"
  else
    ln -sf "${reads[0]}" "$rdir/${lib}.fastq.gz"      # single-end: no _R1 tag (metasnek pairs on _R1)
  fi

  # Drive the real pipeline to the DETECTION target only (skips dedup + polish).
  local target="$work/bsj_results/${lib}.bsj_cor.tsv.gz"
  snakemake \
      -s "$BSJ_FINDER/Snakefile" -d "$BSJ_FINDER" \
      --cores "$THREADS" $SMK_ARGS \
      --config readDirectory="$rdir" genome="$GENOME" outdir="$work" \
      -- "$target"

  # Convert per-read BSJ table -> aggregated BED6, then score.
  python3 "$BENCH/bsj_finder_to_bed.py" --in "$target" \
      --out "$OUT_ABS/$lib.bsj.bed" --merge-tol "$MERGE_TOL"
  python3 "$BENCH/evaluate.py" --tool bed --tool-name bsj_finder --calls "$OUT_ABS/$lib.bsj.bed" \
      --dataset "$DATASET" --library "$lib" --out-prefix "$OUT_ABS/eval_${lib}"
}

run_lib pe150_lowerr  "$DATASET/reads_pe150_lowerr_R1.fastq.gz"  "$DATASET/reads_pe150_lowerr_R2.fastq.gz"
run_lib pe150_higherr "$DATASET/reads_pe150_higherr_R1.fastq.gz" "$DATASET/reads_pe150_higherr_R2.fastq.gz"
run_lib pe100_lowerr  "$DATASET/reads_pe100_lowerr_R1.fastq.gz"  "$DATASET/reads_pe100_lowerr_R2.fastq.gz"
run_lib pe100_higherr "$DATASET/reads_pe100_higherr_R1.fastq.gz" "$DATASET/reads_pe100_higherr_R2.fastq.gz"
run_lib se150_lowerr  "$DATASET/reads_se150_lowerr.fastq.gz"
run_lib se150_higherr "$DATASET/reads_se150_higherr.fastq.gz"
run_lib se100_lowerr  "$DATASET/reads_se100_lowerr.fastq.gz"
run_lib se100_higherr "$DATASET/reads_se100_higherr.fastq.gz"
echo "bsj_finder done -> $OUT"
