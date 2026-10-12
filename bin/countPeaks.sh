#!/bin/bash
#
# countPeaks.sh [--single-end] MAPQ SAF OUT BAM [BAM ...]
#
# Count stage: fragments per merged peak with the vendored featureCounts
# (subread v2.1.1), then Check A (bin/checkCounts.sh) on the result.
# Exits nonzero if either fails, so the run stops at STAGE=Count.
#
#   -p --countReadPairs  count fragments. In 2.x -p alone only says the
#                        input is paired, and every read is counted:
#                        every column doubles, with no error. ChIP-seq
#                        v0.7.0 to v0.9.0 ran 2.0.6 that way.
#   -O                   a fragment over two peaks counts in both
#   -Q MAPQ              10, as before; a no-op, since the postProcess
#                        BAMs are already filtered at MAPQ 30
#   --tmpDir "$TMPDIR"   node-local per-job directory (bin/runStage.sh)
#
#   --single-end         single-end BAMs (pipe.sh -s): count reads, no -p
#                        and no --countReadPairs; Check A in single-end
#                        mode
#
# Defaults that matter: singletons count as fragments (no -B), chimeric
# pairs count (no -C; relevant only to --proper-pair-off runs),
# multimappers do not.
#
# Run from the analysis directory; BAM paths become the column names.
#
# Copied from ATAC-seq devs 28978c9; --single-end is a ChIP-seq addition.
#

set -o pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/loadTools.sh

SINGLE_END=No
if [ "${1:-}" == "--single-end" ]; then
    SINGLE_END=Yes
    shift
fi

if [ "$#" -lt "4" ]; then
    echo
    echo "    usage: countPeaks.sh [--single-end] MAPQ SAF OUT BAM [BAM ...]"
    echo
    exit 1
fi

MAPQ=$1
SAF=$2
OUT=$3
shift 3

THREADS=${SLURM_CPUS_PER_TASK:-1}

$SDIR/featureCounts -v

if [ "$SINGLE_END" == "Yes" ]; then
    PAIR_OPTS=""
    CHECK_OPTS="--single-end"
else
    PAIR_OPTS="-p --countReadPairs"
    CHECK_OPTS=""
fi

$SDIR/featureCounts -O -Q $MAPQ $PAIR_OPTS -T $THREADS \
    --tmpDir "$TMPDIR" \
    -F SAF -a $SAF \
    -o $OUT \
    "$@" \
    || { echo "FATAL ERROR: featureCounts failed"; exit 1; }

#
# The first failing line goes into the FATAL line, so that MESSAGE in the
# run status names the sample (a ChIP-seq change to ATAC's script).
#
if ! CHECK=$($SDIR/checkCounts.sh $CHECK_OPTS $OUT.summary); then
    echo "$CHECK"
    echo "FATAL ERROR: Check A failed: $(echo "$CHECK" | grep -m1 FAIL | tr -s ' ')"
    exit 1
fi
echo "$CHECK"
