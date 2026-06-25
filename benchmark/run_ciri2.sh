#!/usr/bin/env bash
# Run CIRI2 on the tiny BSJ test dataset and evaluate against ground truth.
# Requires: bwa, perl, CIRI2.pl, python3.  See README_benchmark.md.
# Set CIRI2 to the CIRI2.pl path if it is not on PATH:  export CIRI2=/path/to/CIRI2.pl
set -euo pipefail

DATASET="${1:-../dataset}"
OUT="${2:-results/ciri2}"
THREADS="${THREADS:-4}"
GENOME="$DATASET/genome.fa"
GTF="$DATASET/annotation.gtf"
CIRI2="${CIRI2:-CIRI2.pl}"
mkdir -p "$OUT"

command -v bwa >/dev/null 2>&1 || { echo "ERROR: 'bwa' not found"; exit 1; }
command -v perl >/dev/null 2>&1 || { echo "ERROR: 'perl' not found"; exit 1; }
command -v "$CIRI2" >/dev/null 2>&1 || [ -f "$CIRI2" ] || { echo "ERROR: CIRI2.pl not found (set \$CIRI2)"; exit 1; }

# ---- bwa index ----
[ -e "${GENOME}.bwt" ] || bwa index "$GENOME"

run_lib () {
  local lib="$1"; shift
  local reads=("$@")
  echo ">>> CIRI2 :: $lib"
  bwa mem -t "$THREADS" -T 19 "$GENOME" "${reads[@]}" > "$OUT/$lib.sam"
  perl "$CIRI2" -I "$OUT/$lib.sam" -O "$OUT/$lib.ciri" -F "$GENOME" -A "$GTF" -T "$THREADS"
  python3 evaluate.py --tool ciri2 --calls "$OUT/$lib.ciri" \
       --dataset "$DATASET" --library "$lib" --out-prefix "$OUT/eval_${lib}"
  rm -f "$OUT/$lib.sam"   # SAMs are large; drop after parsing
}

run_lib pe_lowerr  "$DATASET/reads_pe_lowerr_R1.fastq.gz"  "$DATASET/reads_pe_lowerr_R2.fastq.gz"
run_lib pe_higherr "$DATASET/reads_pe_higherr_R1.fastq.gz" "$DATASET/reads_pe_higherr_R2.fastq.gz"
run_lib se_lowerr  "$DATASET/reads_se_lowerr.fastq.gz"
run_lib se_higherr "$DATASET/reads_se_higherr.fastq.gz"
echo "CIRI2 done -> $OUT"
