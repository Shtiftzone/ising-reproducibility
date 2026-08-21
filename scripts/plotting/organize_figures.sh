#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

FIGURES_DIR="results/figures"

mkdir -p \
    "$FIGURES_DIR/square/cell" \
    "$FIGURES_DIR/square/vertex" \
    "$FIGURES_DIR/triangular/cell" \
    "$FIGURES_DIR/triangular/vertex"

for COMBO in \
    "square cell" \
    "square vertex" \
    "triangular cell" \
    "triangular vertex"
do
    set -- $COMBO
    LATTICE=$1
    REP=$2

    DEST="$FIGURES_DIR/$LATTICE/$REP"

    files=(
        "$FIGURES_DIR"/*"${LATTICE}_${REP}"*.png
    )

    if [ ${#files[@]} -eq 0 ]; then
        echo "No figures to move for $LATTICE/$REP"
        continue
    fi

    echo "Organizing $LATTICE/$REP..."

    for file in "${files[@]}"; do
        echo "  $(basename "$file")"
        mv "$file" "$DEST/"
    done
done

echo
echo "Figure organization complete."
