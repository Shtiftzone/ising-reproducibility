#!/bin/bash
set -euo pipefail

export LC_ALL=C
export LC_NUMERIC=C

if [[ $# -ne 4 ]]; then
    echo "Usage:"
    echo "  $0 <lattice_type> <representation> <temperatures_dir> <simulation_results_dir>"
    echo
    echo "Examples:"
    echo "  Square, spin as cell:"
    echo "    $0 square cell data/temperatures/square /path/to/square/cell"
    echo
    echo "  Square, spin as vertex:"
    echo "    $0 square vertex data/temperatures/square /path/to/square/vertex"
    echo
    echo "  Triangular, spin as vertex:"
    echo "    $0 triangular vertex data/temperatures/triangular /path/to/triangular/vertex"
    exit 1
fi

LATTICE_TYPE="$1"
REPRESENTATION="$2"
TEMPERATURES_DIR="$3"
RESULTS_DIR="$4"

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

DATA_DIR="data/eul2d/$LATTICE_TYPE/$REPRESENTATION"

mkdir -p "$DATA_DIR"

mapfile -t TEMP_FILES < <(
    find "$TEMPERATURES_DIR" \
        -maxdepth 1 \
        -type f \
        -name 'temperatures_*.txt' \
    | sort -V
)

if [[ ${#TEMP_FILES[@]} -eq 0 ]]; then
    echo "No temperatures_*.txt files found in $TEMPERATURES_DIR"
    exit 1
fi

echo "Lattice type:    $LATTICE_TYPE"
echo "Representation:  $REPRESENTATION"
echo "Results dir:     $RESULTS_DIR"
echo "Output dir:      $DATA_DIR"
echo

echo "Found temperature files:"
printf '  %s\n' "${TEMP_FILES[@]}"
echo

for TEMP_FILE in "${TEMP_FILES[@]}"; do

    SIZE=$(
        basename "$TEMP_FILE" \
        | sed -E 's/^temperatures_([0-9]+)\.txt$/\1/'
    )

    if [[ -z "$SIZE" ]]; then
        echo "Could not read SIZE from filename: $TEMP_FILE"
        continue
    fi

    SIZE_DIR="$DATA_DIR/size_$SIZE"
    mkdir -p "$SIZE_DIR"

    LEFT=$(awk 'NR==1 {print $1}' "$TEMP_FILE")
    RIGHT=$(awk 'END {print $1}' "$TEMP_FILE")

    echo "=============================================="
    echo "Processing L=$SIZE"
    echo "Temperature range: $LEFT - $RIGHT"
    echo "Temperature file: $TEMP_FILE"
    echo "Output directory: $SIZE_DIR"
    echo "=============================================="

    while read -r T IDX; do

        [[ -z "${T:-}" || -z "${IDX:-}" ]] && continue

        TFORMATTED=$(printf "%.4f" "$T")

        SDIR=""

        if [[ -d "$RESULTS_DIR/T_${TFORMATTED}/size_${SIZE}" ]]; then

            SDIR="$RESULTS_DIR/T_${TFORMATTED}/size_${SIZE}"

        else

            MATCH=$(
                find "$RESULTS_DIR" \
                    -maxdepth 1 \
                    -type d \
                    -name "T_${TFORMATTED}*" \
                | sort -V \
                | head -n 1 \
                || true
            )

            if [[ -n "$MATCH" && -d "$MATCH/size_${SIZE}" ]]; then
                SDIR="$MATCH/size_${SIZE}"
            fi
        fi

        if [[ -z "$SDIR" ]]; then
            echo \
                "  Skipping T=$TFORMATTED, idx=$IDX: " \
                "missing directory T_${TFORMATTED}*/size_${SIZE}"
            continue
        fi

        echo "  T=$TFORMATTED  idx=$IDX"
        echo "    source: $SDIR"

        ME="$SDIR/mefile.txt"
        EP="$SDIR/epfile.txt"
        EN="$SDIR/enfile.txt"
        FP="$SDIR/fpfile.txt"
        FN="$SDIR/fnfile.txt"

        MISSING=0

        for f in "$ME" "$EP" "$EN" "$FP" "$FN"; do
            if [[ ! -s "$f" ]]; then
                echo "    missing or empty file: $f"
                MISSING=1
            fi
        done

        [[ $MISSING -eq 1 ]] && continue

        lines_me=$(wc -l < "$ME")
        lines_ep=$(wc -l < "$EP")
        lines_en=$(wc -l < "$EN")
        lines_fp=$(wc -l < "$FP")
        lines_fn=$(wc -l < "$FN")

        MIN_LINES=$lines_me

        (( lines_ep < MIN_LINES )) && MIN_LINES=$lines_ep
        (( lines_en < MIN_LINES )) && MIN_LINES=$lines_en
        (( lines_fp < MIN_LINES )) && MIN_LINES=$lines_fp
        (( lines_fn < MIN_LINES )) && MIN_LINES=$lines_fn

        if (( MIN_LINES == 0 )); then
            echo "    skipping: MIN_LINES = 0"
            continue
        fi

        OUTFILE="$SIZE_DIR/conf-${TFORMATTED}-${IDX}.dat"

        rm -f "$OUTFILE"

        paste \
            <(head -n "$MIN_LINES" "$ME") \
            <(head -n "$MIN_LINES" "$EP") \
            <(head -n "$MIN_LINES" "$EN") \
            <(head -n "$MIN_LINES" "$FP") \
            <(head -n "$MIN_LINES" "$FN") \
        | awk -F'\t' '
            BEGIN {
                OFS=" "
            }
            {
                split($1, me, ",")
                split($4, fp, ",")
                split($5, fn, ",")

                print \
                    me[1], \
                    me[2], \
                    $2, \
                    $3, \
                    fp[1], \
                    fp[2], \
                    fp[3], \
                    fn[1], \
                    fn[2], \
                    fn[3]
            }
        ' > "$OUTFILE"

        echo "    wrote: $OUTFILE  lines: $MIN_LINES"

    done < "$TEMP_FILE"

    echo
done

echo "=============================================="
echo "Done."
echo "=============================================="