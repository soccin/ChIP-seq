#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
PIPENAME="ChIP-Seq"

function usage {
    echo
    echo "usage: $PIPENAME/postPipe.sh ODIR"
    echo "version=$SCRIPT_VERSION"
    echo ""
    echo
    exit
}

if [ "$#" -lt "1" ]; then
    usage
fi

ODIR=$1

mkdir $ODIR/bam
mv $ODIR/*.ba? $ODIR/bam

mkdir $ODIR/bed
mv $ODIR/*.bed.gz $ODIR/bed

Rscript --no-save $SDIR/qc_ChIPSeq_01.R

mkdir -p $ODIR/qc/peaks
mv *___sigPeaks* $ODIR/qc/peaks
mv *___Volcano* $ODIR/qc/peaks
cp qcChIPSeq_* $ODIR/qc/

echo
echo Should create a manifest/group file for stage 2 qc
echo or use _sample_grouping.txt
echo
echo and then run
echo "    Rscript --no-save $SDIR/qc_ChIPSeq_02.R manifest.txt"
echo "    mv *_ChIPSeqQC_*.pdf $ODIR/qc"
echo

if [ -e manifest.txt ]; then
    Rscript --no-save $SDIR/qc_ChIPSeq_02.R manifest.txt
    cp qcChIPSeq_* $ODIR/qc/
    echo
    echo Maybe run homer annotation
    echo    ./ChIP-seq/annotateWithHomer.sh
    echo
fi

