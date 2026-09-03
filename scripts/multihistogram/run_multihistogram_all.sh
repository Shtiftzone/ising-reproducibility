#!/bin/bash
set -euo pipefail

export LC_ALL=C
export LC_NUMERIC=C

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage:"
    echo "  $0 <lattice_type> <representation> <mh_executable_dir> [refinement_factor]"
    echo
    echo "Examples:"
    echo "  Square, cell representation:"
    echo "    $0 square cell src/multihistogram"
    echo
    echo "  Triangular, vertex representation:"
    echo "    $0 triangular vertex src/multihistogram"
    echo
    echo "  With refinement factor 100:"
    echo "    $0 square cell src/multihistogram 100"
    exit 1
fi

LATTICE_TYPE="$1"
REPRESENTATION="$2"
MH_EXECUTABLE_DIR="$3"
REFINEMENT_FACTOR="${4:-50}"

case "$LATTICE_TYPE" in
    square|triangular)
        ;;
    *)
        echo "Invalid lattice type: $LATTICE_TYPE"
        echo "Expected: square or triangular"
        exit 1
        ;;
esac

case "$REPRESENTATION" in
    cell|vertex)
        ;;
    *)
        echo "Invalid representation: $REPRESENTATION"
        echo "Expected: cell or vertex"
        exit 1
        ;;
esac

TEMPERATURES_DIR="data/temperatures/$LATTICE_TYPE"
INPUT_DATA_DIR="data/eul2d/$LATTICE_TYPE/$REPRESENTATION"
OUTPUT_DIR="results/multihistogram/$LATTICE_TYPE/$REPRESENTATION"

RW_EXECUTABLE="$MH_EXECUTABLE_DIR/Rw"

if [[ ! -x "$RW_EXECUTABLE" ]]; then
    echo "Missing or non-executable Rw binary: $RW_EXECUTABLE"
    echo "Build it first, for example:"
    echo "  cd $MH_EXECUTABLE_DIR && make"
    exit 1
fi

if [[ ! -d "$TEMPERATURES_DIR" ]]; then
    echo "Missing temperatures directory: $TEMPERATURES_DIR"
    exit 1
fi

if [[ ! -d "$INPUT_DATA_DIR" ]]; then
    echo "Missing multihistogram input directory: $INPUT_DATA_DIR"
    exit 1
fi

if ! [[ "$REFINEMENT_FACTOR" =~ ^[0-9]+$ ]] || (( REFINEMENT_FACTOR < 1 )); then
    echo "refinement_factor must be a positive integer."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

SIZES=(64 96 128 192 256 384 512 768 1024 1536 2048 3072)

echo "Lattice type:    $LATTICE_TYPE"
echo "Representation:  $REPRESENTATION"
echo "Temperatures:    $TEMPERATURES_DIR"
echo "Input data:      $INPUT_DATA_DIR"
echo "Output:          $OUTPUT_DIR"
echo

for GRID_SIZE in "${SIZES[@]}"; do
    TEMPS_FILE="$TEMPERATURES_DIR/temperatures_${GRID_SIZE}.txt"
    SIZE_DATA_DIR="$INPUT_DATA_DIR/size_${GRID_SIZE}"

    if [[ ! -f "$TEMPS_FILE" ]]; then
        echo "Skipping L=$GRID_SIZE: missing temperature file $TEMPS_FILE"
        continue
    fi

    if [[ ! -d "$SIZE_DATA_DIR" ]]; then
        echo "Skipping L=$GRID_SIZE: missing input directory $SIZE_DATA_DIR"
        continue
    fi

    TMIN=$(awk 'NR==1 {print $1}' "$TEMPS_FILE")
    TMAX=$(awk 'END {print $1}' "$TEMPS_FILE")

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
    echo "Lattice: $LATTICE_TYPE"
    echo "Representation: $REPRESENTATION"
    echo "Temperature file: $TEMPS_FILE"
    echo "Input data directory: $SIZE_DATA_DIR"
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
        "$OUTPUT_FILE" \
        "$INPUT_DATA_DIR"

    echo
done

echo "All available multihistogram runs completed."