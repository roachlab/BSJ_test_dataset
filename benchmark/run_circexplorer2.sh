#!/usr/bin/env bash
# Run CIRCexplorer2 on the tiny BSJ test dataset and evaluate against ground truth.
# Requires: STAR, CIRCexplorer2, gtfToGenePred (UCSC), python3.  See README_benchmark.md.
set -euo pipefail

DATASET="${1:-../dataset}"
OUT="${2:-results/circexplorer2}"
THREADS="${THREADS:-4}"
GENOME="$DATASET/genome.fa"
GTF="$DATASET/annotation.gtf"
mkdir -p "$OUT"

for t in STAR CIRCexplorer2 gtfToGenePred; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: '$t' not found in PATH (see README_benchmark.md)"; exit 1; }
done

# ---- annotation -> CIRCexplorer2 ref (refFlat-like genePred) ----
gtfToGenePred -genePredExt -geneNameAsName2 "$GTF" "$OUT/tmp.genePred"
awk 'BEGIN{OFS="\t"}{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' "$OUT/tmp.genePred" > "$OUT/ref.txt"

# ---- STAR index (SA index size scaled to the small genome) ----
GLEN=$(grep -v '>' "$GENOME" | tr -d '\n' | wc -c)
SAIDX=$(python3 -c "import math;print(min(14,int(math.log2($GLEN)/2-1)))")
if [ ! -e "$OUT/star_index/SAindex" ]; then
  mkdir -p "$OUT/star_index"
  STAR --runMode genomeGenerate --genomeDir "$OUT/star_index" \
       --genomeFastaFiles "$GENOME" \
       --genomeSAindexNbases "$SAIDX" --runThreadN "$THREADS"
fi

run_lib () {
  local lib="$1"; shift
  local overhang="$1"; shift
  local reads=("$@")
  local pfx="$OUT/$lib."
  echo ">>> CIRCexplorer2 :: $lib"
  STAR --runThreadN "$THREADS" --genomeDir "$OUT/star_index" \
       --readFilesIn "${reads[@]}" --readFilesCommand zcat \
       --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimOutType Junctions \
       --chimScoreMin 1 --alignSJoverhangMin 8 \
       --outSAMtype BAM Unsorted --outFileNamePrefix "$pfx"
  CIRCexplorer2 parse -t STAR "${pfx}Chimeric.out.junction" -b "$OUT/$lib.bsj.bed"
  CIRCexplorer2 annotate -r "$OUT/ref.txt" -g "$GENOME" -b "$OUT/$lib.bsj.bed" \
       -o "$OUT/$lib.circ_known.txt" || echo "  (annotate produced no known circles)"
  python3 evaluate.py --tool circexplorer2_bed --calls "$OUT/$lib.bsj.bed" \
       --dataset "$DATASET" --library "$lib" --out-prefix "$OUT/eval_${lib}_denovo"
  if [ -s "$OUT/$lib.circ_known.txt" ]; then
    python3 evaluate.py --tool circexplorer2_known --calls "$OUT/$lib.circ_known.txt" \
       --dataset "$DATASET" --library "$lib" --out-prefix "$OUT/eval_${lib}_known"
  fi
}

run_lib pe150_lowerr  149 "$DATASET/reads_pe150_lowerr_R1.fastq.gz"  "$DATASET/reads_pe150_lowerr_R2.fastq.gz"
run_lib pe150_higherr 149 "$DATASET/reads_pe150_higherr_R1.fastq.gz" "$DATASET/reads_pe150_higherr_R2.fastq.gz"
run_lib pe100_lowerr   99 "$DATASET/reads_pe100_lowerr_R1.fastq.gz"  "$DATASET/reads_pe100_lowerr_R2.fastq.gz"
run_lib pe100_higherr  99 "$DATASET/reads_pe100_higherr_R1.fastq.gz" "$DATASET/reads_pe100_higherr_R2.fastq.gz"
run_lib se150_lowerr  149 "$DATASET/reads_se150_lowerr.fastq.gz"
run_lib se150_higherr 149 "$DATASET/reads_se150_higherr.fastq.gz"
run_lib se100_lowerr   99 "$DATASET/reads_se100_lowerr.fastq.gz"
run_lib se100_higherr  99 "$DATASET/reads_se100_higherr.fastq.gz"
echo "CIRCexplorer2 done -> $OUT"
