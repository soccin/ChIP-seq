#!/bin/bash

# CMD:
#    bsub -n 1 -q control -o LSF.CTRL/ -J CTRL.ChIP ./pipe.sh
#

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
PIPENAME="ChIP-Seq"

##
# Process command args

TAG=q$PIPENAME

COMMAND_LINE=$*

function usage {
    echo
    echo "usage: $PIPENAME/pipe.sh [--proper-pair-off] BAM1 [BAM2 ... BAMN]"
    echo "version=$SCRIPT_VERSION"
    echo ""
    echo
    exit
}

if [ "$#" -lt "1" ]; then
    usage
fi

PROPER_PAIR="Yes"

while :; do
    case $1 in
        --proper-pair-off) PROPER_PAIR="No"
        ;;

        -*)
        echo
        echo "Invalid option ["$1"]"
        usage
        exit
        ;;

        *) break
        ;;

    esac
    shift
done

echo PROPER_PAIR=$PROPER_PAIR
BAMS=$*
echo SDIR=$SDIR
echo BAMS=$BAMS

exit

ODIR=out
mkdir -p $ODIR

RUNTIME="-We 119"
echo $BAMS \
    | xargs -n 1 bsub $RUNTIME -o LSF.POST/ -J ${TAG}_POST2_$$ -R "rusage[mem=24]" \
        echo $SDIR/postMapBamProcessing_ChIPSeq.sh $ODIR $PROPER_PAIR

exit

bSync ${TAG}_POST2_$$

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.BW/ -J ${TAG}_BW2_$$ -R "rusage[mem=24]" \
        $SDIR/makeBigWigFromBEDZ.sh

exit

ls *.bed.gz \
    | xargs -n 1 bsub $RUNTIME -o LSF.CALLP/ -J ${TAG}_CALLP2_$$ -n 3 -R "rusage[mem=24]" \
        $SDIR/callPeaks_ATACSeq.sh

bSync ${TAG}_CALLP2_$$

bsub $RUNTIME -o LSF.CALLP/ -J ${TAG}_MergePeaks_$$ -n 3 -R "rusage[mem=24]" \
    $SDIR/mergePeaksToSAF.sh callpeaks \>macsPeaksMerged.saf

PBAMS=$(ls *_postProcess.bam)
bsub $RUNTIME -o LSF.CALLP/ -J ${TAG}_Count_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_MergePeaks_$$)" \
    $SDIR/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a macsPeaksMerged.saf \
        -o peaks_raw_fcCounts.txt \
        $PBAMS

bsub $RUNTIME -o LSF.DESEQ/ -J ${TAG}_DESEQ_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_Count_$$)" \
    Rscript --no-save $SDIR/getDESeqScaleFactors.R
