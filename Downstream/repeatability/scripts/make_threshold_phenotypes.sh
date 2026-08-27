#!/bin/bash
# make_threshold_phenotypes.sh
# Usage: bash make_threshold_phenotypes.sh <pheno_in> <out> [<low_max>] [<high_min>]
#
# Reads a continuous phenotype file (FID IID ALTITUDE_M) and writes a binary
# phenotype file using hardcoded altitude thresholds:
#   low  = altitude <= LOW_MAX  (default 200m) → 1 (control)
#   high = altitude >= HIGH_MIN (default 1000m) → 2 (case)
#   intermediate                                → -9 (missing/excluded)

set -euo pipefail

PHENO_IN=$1
OUT=$2
LOW_MAX=${3:-200}
HIGH_MIN=${4:-1000}

echo "  Thresholds: low ≤ ${LOW_MAX} m, high ≥ ${HIGH_MIN} m"

awk -v lo="$LOW_MAX" -v hi="$HIGH_MIN" '
    {
        alt = $3
        if      (alt+0 <= lo+0) ph = 1
        else if (alt+0 >= hi+0) ph = 2
        else                    ph = -9
        print $1, $2, ph
    }
' OFS='\t' "$PHENO_IN" > "$OUT"

N_LOW=$(awk '$3 == 1'  "$OUT" | wc -l)
N_HIGH=$(awk '$3 == 2' "$OUT" | wc -l)
N_EXCL=$(awk '$3 == -9' "$OUT" | wc -l)

echo "  Samples: low=${N_LOW}, high=${N_HIGH}, excluded=${N_EXCL}"
echo "  threshold phenotype → $OUT"