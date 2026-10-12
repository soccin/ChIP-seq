#!/bin/bash
#
# runStage.sh CMD [ARG ...]
#
# Run one stage command and append an exit trailer to the job log.
# pipe.sh submits every stage job through this script.
#
# Slurm writes nothing of its own into a job's -o file, so a log alone
# cannot show whether the job worked. The trailer
#
#   #CHIP_EXIT=<rc>
#
# is the last line the job writes. A missing trailer means the job never
# reached the end (killed for memory or time, node failure); then sacct
# is the only record. bin/checkRun.sh reads the trailer as a fallback when
# sacct no longer has the job.
#
# The trailer goes to stderr. The MergePeaks stage redirects the command's
# stdout into macsPeaksMerged.saf; stderr still lands in the job log.
#
# Every temp file of the job goes to a per-job directory on the node's
# /localscratch, never to /tmp (see bin/loadTools.sh). It is removed when
# the job ends, including on scancel or a time limit, when Slurm sends
# SIGTERM to every process of the job before SIGKILL.
#

CHIP_LOCAL_TMP=${CHIP_LOCAL_TMP:-/localscratch/$USER}

mkdir -p "$CHIP_LOCAL_TMP" \
    && TMPDIR=$(mktemp -d -p "$CHIP_LOCAL_TMP" "chip.${SLURM_JOB_ID:-$$}.XXXXXX") \
    || { echo "#CHIP_EXIT=1 cannot create a TMPDIR under $CHIP_LOCAL_TMP" >&2; exit 1; }
export TMPDIR

trap 'rm -rf "$TMPDIR"' EXIT
trap 'exit 143' TERM
trap 'exit 130' INT

echo "#CHIP_HOST=$(hostname) START=$(date '+%Y-%m-%d %H:%M:%S')" >&2
echo "#CHIP_CMD=$*" >&2
echo "#CHIP_TMPDIR=$TMPDIR" >&2

"$@"
RC=$?

echo "#CHIP_EXIT=$RC END=$(date '+%Y-%m-%d %H:%M:%S')" >&2
exit $RC
