#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

GBUILD=$1
BEDZ=$2

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


ODIR=($dirname $BEDZ)/profiles
mkdir -p $ODIR
OUT=$ODIR/$OUTFILE

GENOME=$SDIR/lib/genomes/${GBUILD}.genome

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
    | $SDIR/wigToBigWig stdin $GENOME $OUT

