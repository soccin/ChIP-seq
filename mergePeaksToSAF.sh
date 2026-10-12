#!/bin/bash
#
# mergePeaksToSAF.sh MACS_OUTDIR >SAF
#
# MergePeaks stage: merge every narrowPeak/broadPeak file under
# MACS_OUTDIR (peaks within 500 bp) into one SAF on stdout.
#
# The bedtools merge start is 0-based and is written into the 1-based
# SAF unchanged, so each peak starts one base early (00.ISSUES.md). Kept
# as it is: fixing it changes the counts.
#

set -e -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

if [ "$#" != "1" ]; then
    echo "usage: mergePeaksToSAF.sh MACS_OUTDIR >SAF" >&2
    exit 1
fi

MACS_OUTDIR=$1

PEAKS=$(find $MACS_OUTDIR | egrep "narrowPeak$|broadPeak$") \
    || { echo "FATAL ERROR: no peak files under $MACS_OUTDIR" >&2; exit 1; }

sort -S 16g -T "$TMPDIR" -k1,1V -k2,2n $PEAKS \
    | bedtools merge -i - -d 500 \
    | awk '{print "Peak_"++s,$0,"+"}' | tr ' ' '\t' \
    || { echo "FATAL ERROR: peak merge failed" >&2; exit 1; }
