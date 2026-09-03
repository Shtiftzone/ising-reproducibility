#!/bin/bash
set -euo pipefail
export LC_ALL=C
export LC_NUMERIC=C

if [[ $# -lt 3 || $# -gt 7 ]]; then
    echo "Usage:"
    echo "  $0 <executables_dir> <results_dir> <construction> [T_start] [T_end] [T_step] [Nconf]"
    echo
    echo "Construction:"
    echo "  cell"
    echo "  vertex"
    echo
    echo "Example:"
    echo "  $0 src/simulation results/square_cell_simulations cell 2.22 2.32 0.0005 200"
    echo "  $0 src/simulation results/square_vertex_simulations vertex 2.22 2.32 0.0005 200"
    exit 1
fi

EXECUTABLES_DIR="$1"
RESULTS_DIR="$2"
CONSTRUCTION="$3"
T_START="${4:-2.22}"
T_END="${5:-2.32}"
T_STEP="${6:-0.0005}"
NCONF="${7:-200}"

if [[ "$CONSTRUCTION" != "cell" && "$CONSTRUCTION" != "vertex" ]]; then
    echo "Error: construction must be one of: cell, vertex"
    exit 1
fi

SEED_FILE="${RESULTS_DIR}/seed.txt"
LOG_FILE="${RESULTS_DIR}/simulation_logs_square_${CONSTRUCTION}.txt"

mkdir -p "$RESULTS_DIR"

if [[ ! -f "$SEED_FILE" ]]; then
    echo "123456" > "$SEED_FILE"
fi

command -v awk >/dev/null 2>&1 || { echo "Missing awk."; exit 1; }

TIME_BIN="/usr/bin/time"
if [[ ! -x "$TIME_BIN" ]]; then
    TIME_BIN="$(command -v time || true)"
fi

if [[ -z "${TIME_BIN:-}" ]]; then
    echo "Could not find time command."
    exit 1
fi

SIZES=(64 96 128 192 256 384 512 768 1024 1536 2048 3072)

echo "Starting square-lattice 2D Ising simulations"
echo "Construction: $CONSTRUCTION"
echo "Executables directory: $EXECUTABLES_DIR"
echo "Results directory: $RESULTS_DIR"
echo "Temperature range: $T_START to $T_END with step $T_STEP"
echo "Nconf: $NCONF"
echo "Seed file: $SEED_FILE"
echo

T="$T_START"

while awk -v t="$T" -v end="$T_END" 'BEGIN { exit !(t <= end + 1e-12) }'; do
    TEMP_FMT=$(awk -v t="$T" 'BEGIN { printf "%.5f", t }')
    TEMP_DIR="${RESULTS_DIR}/T_${TEMP_FMT}"
    mkdir -p "$TEMP_DIR"

    echo "=============================================="
    echo "Temperature T=$TEMP_FMT"
    echo "Output directory: $TEMP_DIR"
    echo "=============================================="

    for SIZE in "${SIZES[@]}"; do
        EXEC="${EXECUTABLES_DIR}/ising2d-${CONSTRUCTION}-${SIZE}"

        if [[ ! -x "$EXEC" ]]; then
            echo "Skipping L=$SIZE: missing executable $EXEC"
            continue
        fi

        SIZE_DIR="${TEMP_DIR}/size_${SIZE}"
        mkdir -p "$SIZE_DIR"

        echo "Running square lattice, construction=$CONSTRUCTION, L=$SIZE, T=$TEMP_FMT"

        "$TIME_BIN" -v "$EXEC" "$T" "$NCONF" "$SEED_FILE" \
            "${SIZE_DIR}/mefile.txt" \
            "${SIZE_DIR}/bpfile.txt" \
            "${SIZE_DIR}/bnfile.txt" \
            "${SIZE_DIR}/epfile.txt" \
            "${SIZE_DIR}/enfile.txt" \
            "${SIZE_DIR}/fpfile.txt" \
            "${SIZE_DIR}/fnfile.txt" \
            10 \
            2>> "$LOG_FILE"
    done

    T=$(awk -v t="$T" -v step="$T_STEP" 'BEGIN { printf "%.10f", t + step }')
done

echo "All square-lattice 2D simulations completed."