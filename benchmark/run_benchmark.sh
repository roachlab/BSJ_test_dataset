#!/usr/bin/env bash
# Master benchmark: run CIRCexplorer2 + CIRI2 on all libraries, then summarise.
# Your own pipeline: emit BED6 (chrom,start,end,name,junction_reads,strand) per library and run:
#   python3 evaluate.py --tool bed --calls mypipe_pe150_lowerr.bed --dataset ../dataset \
#           --library pe150_lowerr --out-prefix results/mypipeline/eval_pe150_lowerr
set -euo pipefail
DATASET="${1:-../dataset}"
bash run_circexplorer2.sh "$DATASET" results/circexplorer2 || echo "CIRCexplorer2 step failed/skipped"
bash run_ciri2.sh        "$DATASET" results/ciri2        || echo "CIRI2 step failed/skipped"
bash run_bsj_finder.sh   "$DATASET" results/bsj_finder   || echo "bsj_finder step failed/skipped"
echo; echo "================ BENCHMARK SUMMARY ================"
python3 combine_results.py results
