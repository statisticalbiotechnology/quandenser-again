#!/bin/bash
# Search one set of Quandenser consensus spectra from PXD002370 with Crux Tide
# and rescore with Percolator, the settings every comparison here uses:
# concatenated target-decoy, 20 ppm precursor window, 0.02 Th fragment bins.
#
#   search.sh <results_dir> <out_dir> [tide_index]
#
# The default index is UPS1 + reviewed yeast, built once under
# PXD002370/search; pass another for other datasets.
set -euo pipefail
RES=$(readlink -f "$1"); OUT=$(readlink -f -m "$2")
P=/home/lukask/data/PXD002370
INDEX=$(readlink -f "${3:-$P/search/tide_index}"); IDIR=$(dirname "$INDEX")
IMG=quay.io/biocontainers/crux-toolkit:4.1--h503566f_3
mkdir -p "$OUT"
docker run --rm --memory 12g -v "$IDIR":"$IDIR":ro -v "$RES":"$RES":ro -v "$OUT":"$OUT" -w "$OUT" $IMG \
    crux tide-search --overwrite T --concat T --num-threads 4 \
    --precursor-window 20 --precursor-window-type ppm --mz-bin-width 0.02 --mz-bin-offset 0.0 \
    --output-dir "$OUT/tide" "$RES"/quandenser/consensus_spectra/*.ms2 "$INDEX" > "$OUT/tide.log" 2>&1
docker run --rm --memory 12g -v "$OUT":"$OUT" -w "$OUT" $IMG \
    crux percolator --overwrite T --output-dir "$OUT/percolator" "$OUT/tide/tide-search.txt" > "$OUT/percolator.log" 2>&1
awk -F'\t' 'NR>1 && $8<=0.01 {n++} END {print n " PSMs at 1% FDR"}' "$OUT/percolator/percolator.target.psms.txt"
