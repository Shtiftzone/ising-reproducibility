#!/bin/bash
set -euo pipefail

export LC_ALL=C
export LC_NUMERIC=C

PARAM_FILE="${1:-data/analysis_parameters/reweighting_windows.csv}"
OUTDIR="${2:-data/temperatures/square}"

mkdir -p "$OUTDIR"

if [[ ! -f "$PARAM_FILE" ]]; then
    echo "Missing parameter file: $PARAM_FILE"
    exit 1
fi

echo "Generating square-lattice reweighting temperature files"
echo "Parameter file: $PARAM_FILE"
echo "Output directory: $OUTDIR"

tail -n +2 "$PARAM_FILE" | while IFS=',' read -r LATTICE SIZE TMIN TMAX STEP; do
    [[ "$LATTICE" != "square" ]] && continue

    OUTFILE="$OUTDIR/temperatures_${SIZE}.txt"
    rm -f "$OUTFILE"

    T=$(printf "%.4f" "$TMIN")
    index=0

    while awk -v t="$T" -v r="$TMAX" 'BEGIN { exit !(t <= r + 1e-12) }'; do
        printf "%.4f %d\n" "$T" "$index" >> "$OUTFILE"
        T=$(awk -v t="$T" -v s="$STEP" 'BEGIN { printf "%.4f", t+s }')
        index=$((index + 1))
    done

    echo "L=$SIZE: wrote $OUTFILE with $index temperatures"
done

echo "Done."
