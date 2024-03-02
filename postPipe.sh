#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

source $SDIR/lsf.sh

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

ODIR=$(realpath $1)

if [ ! -e "$ODIR/bam" ]; then
    mkdir $ODIR/bam
    mv $ODIR/*.ba? $ODIR/bam

    mkdir $ODIR/bed
    mv $ODIR/*.bed.gz $ODIR/bed
fi

Rscript --no-save $SDIR/qc_ChIPSeq_01.R

mkdir -p $ODIR/qc/peaks
mv *___sigPeaks* $ODIR/qc/peaks
mv *___Volcano* $ODIR/qc/peaks
cp qcChIPSeq_* $ODIR/qc/

date >01_POST_DONE
cat << EOF | tee -a 01_POST_DONE

Should create a manifest/group file `manifest.txt` for
stage 2 qc or use _sample_grouping.txt

and then run
    Rscript --no-save $SDIR/qc_ChIPSeq_02.R manifest.txt
    mv qcChIPSeq_*.pdf qcChIPSeq_*.xlsx $ODIR/qc

EOF

if [ -e manifest.txt ]; then
    Rscript --no-save $SDIR/qc_ChIPSeq_02.R manifest.txt
    cp qcChIPSeq_* $ODIR/qc/

cat << EOF | tee -a 01_POST_DONE

Maybe run homer annotation (but check if HUMAN or fix!!!)

    find out/macs/Proj_* | egrep "broad|narrow" \
         | xargs -n 1 bsub -o LSF.HOMER/ -J HOMER -W 59 -n 6 \
            ./ChIP-seq/annotateWithHomer.sh

EOF

fi

