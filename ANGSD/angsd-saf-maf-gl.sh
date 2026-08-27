#!/bin/bash

# angsd-saf-maf-gl.sh THREADS REF BAMLIST MININD MINMAF MAXDEPTH DOWNSAMPLE OUTPREFIX
#
# DOWNSAMPLE: fraction to downsample to (e.g. 0.5), or "none" to skip -downSample entirely

set -euo pipefail
 
THREADS=$1
REF=$2
BAMLIST=$3
MININD=$4
MINMAF=$5
MAXDEPTH=$6
DOWNSAMPLE=$7
GENOMINDEPTH=$8
OUTPREFIX=$9
 
DOWNSAMPLE_ARGS=()
if [[ "${DOWNSAMPLE}" != "none" ]]; then
    DOWNSAMPLE_ARGS=(-downSample "${DOWNSAMPLE}")
fi
 
GENO_MIN_DEPTH_ARGS=()
if [[ "${GENOMINDEPTH}" != "none" ]]; then
    GENO_MIN_DEPTH_ARGS=(-geno_minDepth "${GENOMINDEPTH}")
fi
 
echo "N individuals in bam list: $(wc -l < "${BAMLIST}")"
echo "minInd: ${MININD}"
echo "minMaf: ${MINMAF}"
echo "maxDepth / setMaxDepth: ${MAXDEPTH}"
echo "downSample: ${DOWNSAMPLE}"
echo "geno_minDepth: ${GENOMINDEPTH}"
 
angsd -P "${THREADS}" \
    -doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
    -doGeno 2 -doPost 1 -doDepth 1 -maxDepth "${MAXDEPTH}" -dumpCounts 2 \
    -anc "${REF}" -remove_bads 1 -minMapQ 30 -minQ 20 \
    -minInd "${MININD}" -minMaf "${MINMAF}" -setMaxDepth "${MAXDEPTH}" \
    "${DOWNSAMPLE_ARGS[@]}" \
    "${GENO_MIN_DEPTH_ARGS[@]}" \
    -b "${BAMLIST}" \
    -out "${OUTPREFIX}"
 
echo "Done."