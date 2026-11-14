#!/bin/bash
# Synthetic end-to-end test: generate data, run the Xenium pipeline on the Xenium half, this
# pipeline on the CosMx half, then the comparison (10). Nothing is written into the repositories:
# both code folders are copied next to the data and get the synthetic config there.
#
#   tests/run_synthetic.sh <out_dir> <path to the Xenium repository>/experiment1_tma
#   SYN_CROWD_K=Inf tests/run_synthetic.sh ...     # negative control: no crowding simulated
#
# Takes a few minutes. Look at <out_dir>/cosmx_project/tables/10_summary.csv and figures/10_*.png.
set -euo pipefail
OUT=$(mkdir -p "$1" && cd "$1" && pwd)
XEN_CODE=$(cd "$2" && pwd)
HERE=$(cd "$(dirname "$0")/.." && pwd)

Rscript "$HERE/tests/make_synthetic_data.R" "$OUT"

rm -rf "$OUT/code_xenium" "$OUT/code_cosmx"
cp -r "$XEN_CODE" "$OUT/code_xenium" && cp "$OUT"/config_xenium/*.csv "$OUT/code_xenium/config/"
# the Xenium repository's config.R names the clusters of ITS data set; here there are none yet
printf '\nCLUSTER_TO_CELLTYPE <- NULL\nCLUSTER_TO_LINEAGE <- NULL\nINSPECT_CELLPOSE <- FALSE\nINSPECT_N_PER_SLIDE <- 0\n' >> "$OUT/code_xenium/config/config.R"
(cd "$OUT/code_xenium" && export EXP1_PROJECT_DIR="$OUT/xenium_project" &&
   Rscript 01_xenium_object_prep.R && Rscript 02_tma_dearray.R && Rscript 04_annotation_probe_based.R)

cp -r "$HERE" "$OUT/code_cosmx" && cp "$OUT"/config_cosmx/*.csv "$OUT/code_cosmx/config/"
cd "$OUT/code_cosmx"
export EXP1_PROJECT_DIR="$OUT/cosmx_project" EXP1_XENIUM_PROJECT_DIR="$OUT/xenium_project"
for s in 00_cosmx_inspect_input.R 01_cosmx_object_prep.R 02_tma_dearray.R 03_segmentation_qc.R 04_annotation_wtx.R \
         10_xenium_vs_cosmx.R 05_tma_first_insights.R; do
  echo "=== $s"; Rscript "$s"
done
echo "Done. Summary: $OUT/cosmx_project/tables/10_summary.csv"
