#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

export PYTHONPATH=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/lib/python2.7/site-packages/:$PYTHONPATH
MACS=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/bin/macs2

GBUILD=$1
BEDZ=$2

GENOME=$SDIR/lib/genomes/${GBUILD}.genome
if [ ! -e "$GENOME" ]; then

    echo
    echo "  FATAL ERROR: MISSING GENOME OUTFILE"
    echo
    echo "  GBUILD=$GBUILD, GENOME=$GENOME"
    echo

    exit 1

fi

case $GBUILD in
    mm10*)
    MACS_GENOME=mm
    ;;

    b37)
    MACS_GENOME=hs
    ;;

    *)
    echo "Invalid Genome [$GBUILD]"
    exit -1
    ;;
esac

echo MACS_GENOME=$MACS_GENOME

if [ "$#" == "3" ]; then
    scaleFactor=$3
    echo "$BEDZ sizeFactorNorm scaleFactor "$scaleFactor
    OUTFILE=$(basename $BEDZ | sed 's/.bed.gz/.sizeFactorNorm.bw/')
else
    count=$(zcat $BEDZ | cut -f4 | sort -S20g | uniq | wc -l)
    scaleFactor=$(bc -l <<< "10000000/$count")
    echo "$BEDZ 10mNorm scaleFactor "$scaleFactor
    OUTFILE=$(basename $BEDZ | sed 's/.bed.gz/.10mNorm.bw/')
fi

OUTPREDICT=$(basename $BEDZ | sed 's/.bed.gz//')

ODIR=$(dirname $BEDZ)/profiles
mkdir -p $ODIR
OUT=$ODIR/$OUTFILE

# TDIR=/scratch/socci
# mkdir -p $TDIR
# TMP=$(mktemp -p $TDIR)

#
# N.B. -scale argument is multiplicative
#
#    scaledCounts = rawCounts * $scaleFactor
#

zcat $BEDZ \
    | bedtools slop -i - -g $GENOME -s -l 0 -r 0 \
    | bedtools genomecov -i - -g $GENOME -bg -scale $scaleFactor \
    | $SDIR/bin/wigToBigWig stdin $GENOME $OUT

mkdir -p $ODIR/$OUTPREDICT
$MACS predictd -g $MACS_GENOME -i $BEDZ --outdir $ODIR/$OUTPREDICT 2> $ODIR/${OUTPREDICT}.log
