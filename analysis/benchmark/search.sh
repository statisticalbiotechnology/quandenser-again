#!/bin/bash
# Search one set of Quandenser consensus spectra with Crux Tide and rescore
# with Percolator, the settings every comparison here uses: concatenated
# target-decoy, 20 ppm precursor window, 0.02 Th fragment bins.
#
#   search.sh <results_dir> <out_dir> [tide_index]
#
# The default index is UPS1 + reviewed yeast, built once under
# PXD002370/search; pass another for other datasets.
#
# Tide 4.1 does not release memory while it searches: on a file of 200,000
# PASEF consensus spectra (PXD028735) it passed 16 GB. Each consensus file is
# therefore searched in a call of its own, and the outputs are concatenated
# (scan numbers are unique across Quandenser's consensus files, and every row
# names its file). SEARCH_MEMORY sets the container limit (default 24g).
set -euo pipefail
RES=$(readlink -f "$1"); OUT=$(readlink -f -m "$2")
P=/home/lukask/data/PXD002370
INDEX=$(readlink -f "${3:-$P/search/tide_index}"); IDIR=$(dirname "$INDEX")
IMG=quay.io/biocontainers/crux-toolkit:4.1--h503566f_3
MEM=${SEARCH_MEMORY:-24g}
mkdir -p "$OUT/tide"
: > "$OUT/tide.log"
first=1
for MS2 in "$RES"/quandenser/consensus_spectra/*.ms2; do
  PART="$OUT/tide/$(basename "$MS2" .ms2)"
  docker run --rm --memory $MEM -v "$IDIR":"$IDIR":ro -v "$RES":"$RES":ro -v "$OUT":"$OUT" -w "$OUT" $IMG \
      crux tide-search --overwrite T --concat T --num-threads 4 \
      --precursor-window 20 --precursor-window-type ppm --mz-bin-width 0.02 --mz-bin-offset 0.0 \
      --output-dir "$PART" "$MS2" "$INDEX" >> "$OUT/tide.log" 2>&1
  if [ $first -eq 1 ]; then cat "$PART/tide-search.txt"; first=0; else tail -n +2 "$PART/tide-search.txt"; fi
  rm "$PART/tide-search.txt"
done > "$OUT/tide/tide-search.txt"
docker run --rm --memory $MEM -v "$OUT":"$OUT" -w "$OUT" $IMG \
    crux percolator --overwrite T --output-dir "$OUT/percolator" "$OUT/tide/tide-search.txt" > "$OUT/percolator.log" 2>&1
awk -F'\t' 'NR>1 && $8<=0.01 {n++} END {print n " PSMs at 1% FDR"}' "$OUT/percolator/percolator.target.psms.txt"
