#!/bin/bash
#
# makeBigWigFromBEDZ.sh GENOME BEDZ [SCALE_FACTOR]
#
# BW stage. From a clean.bed.gz written by POST, makes under
# <dir of BEDZ>/profiles/:
#
#   <name>.10mNorm.bw          scaled to 10M reads (default), or
#   <name>.sizeFactorNorm.bw   scaled by SCALE_FACTOR when it is given
#   <name>.log                 macs2 predictd output; the fragment length
#                              estimate pipe.sh reads
#
# A predictd failure, or a log with no estimate, fails the job (spec
# decision D11): the run stops at STAGE=BW with this script's FATAL line
# in MESSAGE.
#

set -e -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

. $SDIR/venv/bin/activate
MACS=$SDIR/venv/bin/macs2

if [ "$#" != "2" ] && [ "$#" != "3" ]; then
    echo
    echo "    usage: makeBigWigFromBEDZ.sh GENOME BEDZ [SCALE_FACTOR]"
    echo
    exit 1
fi

GBUILD=$1
BEDZ=$2

GENOME=$($SDIR/bin/genomeInfo.sh $GBUILD GENOME_FILE) \
    || { echo "FATAL ERROR: no genome file for [$GBUILD]"; exit 1; }
MACS_GENOME=$($SDIR/bin/genomeInfo.sh $GBUILD MACS_GSIZE) \
    || { echo "FATAL ERROR: no MACS genome size for [$GBUILD]"; exit 1; }

echo GBUILD=$GBUILD GENOME=$GENOME MACS_GENOME=$MACS_GENOME

if [ "$#" == "3" ]; then
    scaleFactor=$3
    echo "$BEDZ sizeFactorNorm scaleFactor "$scaleFactor
    OUTFILE=$(basename $BEDZ .bed.gz).sizeFactorNorm.bw
else
    count=$(zcat $BEDZ | cut -f4 | sort -S20g -T "$TMPDIR" | uniq | wc -l) \
        || { echo "FATAL ERROR: read count of $BEDZ failed"; exit 1; }
    if [ "$count" -le 0 ]; then
        echo "FATAL ERROR: no reads in $BEDZ"
        exit 1
    fi
    scaleFactor=$(bc -l <<< "10000000/$count")
    echo "$BEDZ 10mNorm count $count scaleFactor "$scaleFactor
    OUTFILE=$(basename $BEDZ .bed.gz).10mNorm.bw
fi

OUTPREDICT=$(basename $BEDZ .bed.gz)

ODIR=$(dirname $BEDZ)/profiles
mkdir -p $ODIR
OUT=$ODIR/$OUTFILE

#
# N.B. -scale argument is multiplicative
#
#    scaledCounts = rawCounts * $scaleFactor
#

zcat $BEDZ \
    | bedtools slop -i - -g $GENOME -s -l 0 -r 0 \
    | bedtools genomecov -i - -g $GENOME -bg -scale $scaleFactor \
    | $SDIR/bin/wigToBigWig stdin $GENOME $OUT \
    || { echo "FATAL ERROR: bigWig build failed for $BEDZ"; exit 1; }

mkdir -p $ODIR/$OUTPREDICT
$MACS predictd \
    -g $MACS_GENOME -i $BEDZ --outdir $ODIR/$OUTPREDICT \
    >$ODIR/${OUTPREDICT}.log 2>&1 \
    || { cat $ODIR/${OUTPREDICT}.log; echo "FATAL ERROR: predictd failed for $BEDZ"; exit 1; }

cat $ODIR/${OUTPREDICT}.log

#
# predictd exits 0 when it cannot build its model ("Too few paired
# peaks"), so the exit status is not enough: no estimate in the log is a
# failure too.
#
if ! grep -q -E "# predicted fragment length is [0-9]+ bps" $ODIR/${OUTPREDICT}.log; then
    echo "FATAL ERROR: predictd gave no fragment length estimate for $BEDZ"
    exit 1
fi
