#!/bin/bash
# make_categorical_phenotypes.sh
# Usage: bash make_categorical_phenotypes.sh <pheno_in> <b2_out>
#
# Reads a continuous phenotype file (FID IID ALTITUDE_M) and writes a categorical
# phenotype file where:
#   low  = bottom 25% of altitude values for this species → 1
#   high = top 25% of altitude values for this species    → 2 
#   intermediate (middle 50%)                             → -9 (excluded)
#
# Thresholds are computed per-species from the input file itself.

set -euo pipefail

PHENO_IN=$1
B2_OUT=$2

# Compute per-species 25th and 75th percentile thresholds using R
# (awk alone doesn't do quantiles cleanly)
THRESHOLDS=$(Rscript - "$PHENO_IN" << 'REOF'
args    <- commandArgs(trailingOnly = TRUE)
pheno   <- read.table(args[1], header = FALSE, col.names = c("FID","IID","ALT"))
q25     <- quantile(pheno$ALT, 0.25, na.rm = TRUE)
q75     <- quantile(pheno$ALT, 0.75, na.rm = TRUE)
cat(q25, q75, "\n")
REOF
)

LOW_MAX=$(echo "$THRESHOLDS" | awk '{print $1}')
HIGH_MIN=$(echo "$THRESHOLDS" | awk '{print $2}')

echo "  Thresholds for $(basename $PHENO_IN): low ≤ ${LOW_MAX} m, high ≥ ${HIGH_MIN} m"

# Write categorical phenotype: 1=low, 2=high, -9=excluded
awk -v lo="$LOW_MAX" -v hi="$HIGH_MIN" '
    {
        alt = $3
        if      (alt <= lo) ph = 1
        else if (alt >= hi) ph = 2
        else                ph = -9
        print $1, $2, ph
    }
' OFS='\t' "$PHENO_IN" > "$B2_OUT"

# Report sample counts
N_LOW=$(awk '$3 == 1' "$B2_OUT" | wc -l)
N_HIGH=$(awk '$3 == 2' "$B2_OUT" | wc -l)
N_EXCL=$(awk '$3 == -9' "$B2_OUT" | wc -l)
echo "  Samples: low=${N_LOW}, high=${N_HIGH}, excluded=${N_EXCL}"
echo "  categorical2 → $B2_OUT"