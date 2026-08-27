#!/bin/bash

# make-sites-list.sh MAFS_GZ OUT_SITES
#
# Extracts chromo, position, major, minor from an ANGSD .mafs.gz file
# (already filtered by -minMaf/-minInd at the angsd step) and builds
# an indexed angsd sites file.

set -euo pipefail

MAFS_GZ=$1
OUT_SITES=$2

zcat "${MAFS_GZ}" \
    | tail -n +2 \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
    > "${OUT_SITES}"

angsd sites index "${OUT_SITES}"

echo "Done."