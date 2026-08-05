#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LC_NUMERIC=C

OUTDIR="${1:-data/temperatures}"
mkdir -p "$OUTDIR"

# Temperature windows used for the multihistogram reweighting analysis.
# These ranges were selected manually around the relevant peak region.
declare -A LEFT
declare -A RIGHT

LEFT[64]=2.2200      ; RIGHT[64]=2.2600
LEFT[96]=2.2275      ; RIGHT[96]=2.2700
LEFT[128]=2.2350     ; RIGHT[128]=2.2800
LEFT[192]=2.2400     ; RIGHT[192]=2.2800
LEFT[256]=2.2450     ; RIGHT[256]=2.2750
LEFT[384]=2.2500     ; RIGHT[384]=2.2750
LEFT[512]=2.2550     ; RIGHT[512]=2.2750
LEFT[768]=2.2600     ; RIGHT[768]=2.2750
LEFT[1024]=2.2625    ; RIGHT[1024]=2.2725
LEFT[1536]=2.2650    ; RIGHT[1536]=2.2725
LEFT[2048]=2.2660    ; RIGHT[2048]=2.2720
LEFT[3072]=2.2670    ; RIGHT[3072]=2.2710

SIZES=(64 96 128 192 256 384 512 768 1024 1536 2048 3072)

STEP=0.0005

echo "Generating reweighting temperature files in: $OUTDIR"

for SIZE in "${SIZES[@]}"; do
    OUTFILE="$OUTDIR/temperatures_${SIZE}.txt"
    L=${LEFT[$SIZE]}
    R=${RIGHT[$SIZE]}

    rm -f "$OUTFILE"

    T=$(printf "%.4f" "$L")
    index=0

    while awk -v t="$T" -v r="$R" 'BEGIN{exit !(t <= r + 1e-12)}'; do
        printf "%.4f %d\n" "$T" "$index" >> "$OUTFILE"
        T=$(awk -v t="$T" -v s="$STEP" 'BEGIN{printf "%.4f", t+s}')
        index=$((index + 1))
    done

    echo "L=$SIZE: wrote $OUTFILE with $index temperatures"
done

echo "Done."
