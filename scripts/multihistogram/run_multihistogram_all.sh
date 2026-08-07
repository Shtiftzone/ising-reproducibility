#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LC_NUMERIC=C

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage:"
    echo "  $0 <mh_executable_dir> <temperatures_dir> <output_dir> [refinement_factor]"
    echo
    echo "Example:"
    echo "  $0 src/multihistogram data/temperatures results/multihistogram"
    echo "  $0 src/multihistogram data/temperatures results/multihistogram 50"
    exit 1
fi

MH_EXECUTABLE_DIR="$1"
TEMPERATURES_DIR="$2"
OUTPUT_DIR="$3"
REFINEMENT_FACTOR="${4:-50}"

RW_EXECUTABLE="$MH_EXECUTABLE_DIR/Rw"

if [[ ! -x "$RW_EXECUTABLE" ]]; then
    echo "Missing or non-executable Rw binary: $RW_EXECUTABLE"
    echo "Build it first, for example:"
    echo "  cd $MH_EXECUTABLE_DIR && make"
    exit 1
fi

if ! [[ "$REFINEMENT_FACTOR" =~ ^[0-9]+$ ]] || (( REFINEMENT_FACTOR < 1 )); then
    echo "refinement_factor must be a positive integer."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SIZES=(64 96 128 192 256 384 512 768 1024 1536 2048 3072)

for GRID_SIZE in "${SIZES[@]}"; do
    TEMPS_FILE="$TEMPERATURES_DIR/temperatures_${GRID_SIZE}.txt"

    if [[ ! -f "$TEMPS_FILE" ]]; then
        echo "Skipping L=$GRID_SIZE: missing temperature file $TEMPS_FILE"
        continue
    fi

    TMIN=$(awk 'NR==1{print $1}' "$TEMPS_FILE")
    TMAX=$(awk 'END{print $1}' "$TEMPS_FILE")

    N_TEMPS=$(wc -l < "$TEMPS_FILE")
    if (( N_TEMPS < 2 )); then
        echo "Skipping L=$GRID_SIZE: too few temperatures in $TEMPS_FILE"
        continue
    fi

    POINTS=$(( (N_TEMPS - 1) * REFINEMENT_FACTOR + 1 ))

    RESULTS_FILE="$OUTPUT_DIR/results_${GRID_SIZE}.dat"
    OUTPUT_FILE="$OUTPUT_DIR/output2_${GRID_SIZE}.dat"

    echo "=============================================="
    echo "Running multihistogram analysis for L=$GRID_SIZE"
    echo "Temperature file: $TEMPS_FILE"
    echo "Tmin=$TMIN"
    echo "Tmax=$TMAX"
    echo "N_temperatures=$N_TEMPS"
    echo "Refinement factor=$REFINEMENT_FACTOR"
    echo "Reweighting points=$POINTS"
    echo "Results file: $RESULTS_FILE"
    echo "Output file: $OUTPUT_FILE"
    echo "=============================================="

    ulimit -s unlimited || true

    "$RW_EXECUTABLE" \
        "$TEMPS_FILE" \
        "$GRID_SIZE" \
        2 \
        "$RESULTS_FILE" \
        "$TMIN" \
        "$TMAX" \
        "$POINTS" \
        "$OUTPUT_FILE"

    echo
done

echo "All available multihistogram runs completed."
