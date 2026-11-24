#!/bin/bash
# Search the raw MS2 chunks of one run, as written by pasef_extract_ms2.py
# with the instrument's precursor and charge, with the settings of search.sh,
# then Percolator. Tide needs about 11 GB for a chunk of 60,000 PASEF spectra.
#   search_raw.sh <chunk_dir> <out_dir> <tide_index> [extra tide-search args]
set -euo pipefail
export SUPPRESS_BOLTDB_WARNING=1
IN=$(readlink -f "$1"); OUT=$(readlink -f -m "$2"); INDEX=$(readlink -f "$3"); shift 3
IDIR=$(dirname "$INDEX")
IMG=quay.io/biocontainers/crux-toolkit:4.1--h503566f_3
mkdir -p "$OUT/tide"; : > "$OUT/tide.log"
first=1
for MS2 in "$IN"/raw.part*.ms2; do
  PART="$OUT/tide/$(basename "$MS2" .ms2)"
  docker run --rm --memory 24g -v "$IDIR":"$IDIR":ro -v "$IN":"$IN":ro -v "$OUT":"$OUT" -w "$OUT" $IMG \
      crux tide-search --overwrite T --concat T --num-threads 4 \
      --precursor-window 20 --precursor-window-type ppm --mz-bin-width 0.02 --mz-bin-offset 0.0 "$@" \
      --output-dir "$PART" "$MS2" "$INDEX" >> "$OUT/tide.log" 2>&1
  if [ $first -eq 1 ]; then cat "$PART/tide-search.txt"; first=0; else tail -n +2 "$PART/tide-search.txt"; fi
  rm "$PART/tide-search.txt"
done > "$OUT/tide/tide-search.txt"
docker run --rm --memory 24g -v "$OUT":"$OUT" -w "$OUT" $IMG \
    crux percolator --overwrite T --output-dir "$OUT/percolator" "$OUT/tide/tide-search.txt" > "$OUT/percolator.log" 2>&1
awk -F'\t' 'NR>1 && $8<=0.01 {n++} END {print n " PSMs at 1% FDR"}' "$OUT/percolator/percolator.target.psms.txt"
echo RAWSEARCH_DONE
