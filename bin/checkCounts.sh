#!/bin/bash
#
# checkCounts.sh [--single-end] SUMMARY [ASSIGNED ...]
#
# Check a featureCounts run before its counts are used. Exits 0 if every
# check passes, 1 if any fails. Prints one line per sample and check,
# whether it passes or not, so the log records the numbers.
#
# Run it from the directory featureCounts ran in (the analysis
# directory): the BAM paths in the SUMMARY header are relative to it.
#
# Check A (pairing), always. For each BAM column of SUMMARY
# (peaks_raw_fcCounts.txt.summary):
#
#   frag   the sum of every row of the column: each fragment is counted
#          once in some status row, assigned or not. Rows are summed, not
#          picked by name; featureCounts 2.x renamed and added rows.
#   pairs  mapped reads / 2, from samtools idxstats (needs the .bai)
#
# and frag / pairs must be within CHECKA_MIN to CHECKA_MAX. A correct
# count gives 1.000 to 1.004 (the excess is singletons whose mate failed
# the MAPQ filter). A BAM whose reads were not paired gives 2.0: the
# featureCounts 1.6.0 bug, or 2.x run with -p but without
# --countReadPairs.
#
# Check B (independent count), when ASSIGNED files are given. Each is the
# one-line output of bin/assignedCheck.sh (BAM<TAB>FRAGMENTS), a count of
# the fragments overlapping the atlas made with samtools alone. It must
# equal the Assigned row of the BAM's column exactly; on the first
# full-size project (12 b38 samples) and the 11 downsampled BAMs they
# were equal for every sample. Every BAM column must have one.
#
# ChIP-seq additions to ATAC-seq devs 28978c9:
#
#   --single-end  for single-end BAMs (pipe.sh -s), counted as reads:
#                 the ratio is frag / mapped reads, same band.
#
#   out/bam/      postPipe.sh moves the BAMs from out/ to out/bam/. If a
#                 BAM in the SUMMARY header does not exist and
#                 <dir>/bam/<name> does, that file is used for Check A,
#                 and the line says so. Check B results are still
#                 matched by the path in the header.
#

set -o pipefail

SDIR=$(cd "$(dirname "$0")" && pwd)

. $SDIR/loadTools.sh

CHECKA_MIN=0.95
CHECKA_MAX=1.05

SINGLE_END=No
if [ "${1:-}" == "--single-end" ]; then
    SINGLE_END=Yes
    shift
fi

if [ "$#" -lt "1" ]; then
    echo
    echo "    usage: checkCounts.sh [--single-end] SUMMARY [ASSIGNED ...]"
    echo
    exit 1
fi

#
# Paired-end: one fragment per two mapped reads. Single-end: one per read.
#
if [ "$SINGLE_END" == "Yes" ]; then
    READS_PER_FRAG=1
    UNIT=reads
else
    READS_PER_FRAG=2
    UNIT=pairs
fi

SUMMARY=$1
shift

if [ ! -s "$SUMMARY" ]; then
    echo "checkCounts: FAILED no summary file [$SUMMARY]"
    exit 1
fi

#
# Check B results by BAM, keyed by absolute path so that the same BAM
# named two ways still matches.
#
declare -A ASSIGNED_BT
for F in "$@"; do
    if [ ! -s "$F" ]; then
        echo "checkCounts: FAILED empty or missing Check B result [$F]"
        exit 1
    fi
    while IFS=$'\t' read -r BAM N; do
        ASSIGNED_BT[$(readlink -m "$BAM")]=$N
    done <"$F"
done

NCOL=$(head -1 "$SUMMARY" | awk -F'\t' '{print NF}')
NBAM=$((NCOL - 1))
NFAIL=0

for ((i = 2; i <= NCOL; i++)); do
    BAM=$(head -1 "$SUMMARY" | cut -f$i)
    SID=$(basename "$BAM" .bam)
    SID=${SID%_postProcess}

    FRAG=$(awk -F'\t' -v c=$i 'NR > 1 {s += $c} END {printf "%d", s}' "$SUMMARY")
    ASSIGNED=$(awk -F'\t' -v c=$i '$1 == "Assigned" {print $c}' "$SUMMARY")

    #
    # Check A
    #
    IBAM=$BAM
    MOVED=""
    if [ ! -e "$BAM" ] && [ -e "$(dirname "$BAM")/bam/$(basename "$BAM")" ]; then
        IBAM=$(dirname "$BAM")/bam/$(basename "$BAM")
        MOVED="  (read from $IBAM)"
    fi

    if READS=$(samtools idxstats "$IBAM" | awk '{s += $3} END {printf "%d", s}') \
           && [ "${READS:-0}" -gt 0 ]; then
        if ! LINE=$(awk -v sid="$SID" -v f="$FRAG" -v r="$READS" -v d=$READS_PER_FRAG \
                        -v unit=$UNIT -v lo=$CHECKA_MIN -v hi=$CHECKA_MAX 'BEGIN {
                ratio = f / (r / d)
                ok = (ratio >= lo && ratio <= hi)
                printf "CHECK_A  %-20s frag=%d %s=%.1f ratio=%.4f  %s\n", \
                    sid, f, unit, r / d, ratio, (ok ? "OK" : "FAIL (outside " lo " to " hi ")")
                exit !ok
            }'); then
            NFAIL=$((NFAIL + 1))
        fi
    else
        LINE="CHECK_A  $(printf '%-20s' $SID) FAIL (no mapped reads from samtools idxstats $IBAM)"
        NFAIL=$((NFAIL + 1))
    fi
    echo "$LINE$MOVED"

    #
    # Check B
    #
    if [ "$#" -gt "0" ]; then
        BT=${ASSIGNED_BT[$(readlink -m "$BAM")]}
        if [ -z "$BT" ]; then
            echo "CHECK_B  $(printf '%-20s' $SID) FAIL (no Check B result for $BAM)"
            NFAIL=$((NFAIL + 1))
        elif [ "$BT" == "$ASSIGNED" ]; then
            echo "CHECK_B  $(printf '%-20s' $SID) assigned=$ASSIGNED samtools=$BT diff=0  OK"
        else
            echo "CHECK_B  $(printf '%-20s' $SID) assigned=$ASSIGNED samtools=$BT diff=$((ASSIGNED - BT))  FAIL (must be equal)"
            NFAIL=$((NFAIL + 1))
        fi
    fi
done

#
# A Check B result for a BAM that is not in SUMMARY means a sample is
# missing from the count matrix.
#
if [ "$#" -gt "0" ] && [ "${#ASSIGNED_BT[@]}" != "$NBAM" ]; then
    echo "CHECK_B  FAIL (${#ASSIGNED_BT[@]} Check B results for $NBAM BAM column(s) in $SUMMARY)"
    NFAIL=$((NFAIL + 1))
fi

if [ "$NFAIL" -gt "0" ]; then
    echo "checkCounts: FAILED $NFAIL check(s), $NBAM BAM(s), $SUMMARY"
    exit 1
fi

echo "checkCounts: OK $NBAM BAM(s), $SUMMARY"
