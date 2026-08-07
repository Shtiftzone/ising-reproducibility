#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LC_NUMERIC=C

OUTDIR="${1:-data/temperatures/triangular}"
mkdir -p "$OUTDIR"

# Temperature windows used for the triangular-lattice
# multihistogram reweighting analysis.
#
# These ranges were selected manually around the relevant peak region.

declare -A LEFT
declare -A RIGHT

LEFT[64]=3.6800      ; RIGHT[64]=3.7200
LEFT[96]=3.6650      ; RIGHT[96]=3.6900
LEFT[128]=3.6550     ; RIGHT[128]=3.6800
LEFT[192]=3.6500     ; RIGHT[192]=3.6700
LEFT[256]=3.6400     ; RIGHT[256]=3.6800
LEFT[384]=3.6400     ; RIGHT[384]=3.6700
LEFT[512]=3.6400     ; RIGHT[512]=3.6600
LEFT[768]=3.6400     ; RIGHT[768]=3.6600
LEFT[1024]=3.6400    ; RIGHT[1024]=3.6550
LEFT[1536]=3.6350    ; RIGHT[1536]=3.6500
LEFT[2048]=3.6350    ; RIGHT[2048]=3.6500
LEFT[3072]=3.6350    ; RIGHT[3072]=3.6500

SIZES=(64 96 128 192 256 384 512 768 1024 1536 2048 3072)

STEP=0.0005

echo "Generating triangular-lattice reweighting temperature files in: $OUTDIR"

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
