#!/bin/bash
# salmon_quant.sh
# Usage: bash salmon_quant.sh <sample> <r1> <r2> <index> <quants_dir> <threads>

set -euo pipefail

SAMPLE=$1
R1=$2
R2=$3
INDEX=$4
QUANTS_DIR=$5
THREADS=$6

echo "Sample   : $SAMPLE"
echo "Start    : $(date)"

mkdir -p "$QUANTS_DIR"

salmon quant \
    -i "$INDEX" \
    -l A \
    -1 "$R1" \
    -2 "$R2" \
    -p "$THREADS" \
    --validateMappings \
    -o "${QUANTS_DIR}/${SAMPLE}"

echo "End: $(date)"