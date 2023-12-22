#!/bin/bash

# CMD:
#    bsub -n 1 -q control -o LSF.CTRL/ -J CTRL.ChIP ./pipe.sh
#

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long)
PIPENAME="ChIP-Seq"

module load bedtools

##
# Process command args

TAG=q$PIPENAME

COMMAND_LINE=$*

function usage {
    echo
    echo "usage: $PIPENAME/pipe.sh [-n|--narrow-peaks] [-o|--outdir <DIR>] [--proper-pair-off] [-s|--single-end-on] --pairing-file <PAIRS> BAM1 [BAM2 ... BAMN]"
    echo "version=$SCRIPT_VERSION"
    echo ""
    echo
    exit
}

if [ "$#" -lt "1" ]; then
    usage
fi

PROPER_PAIR="Yes"
SE="No"
ODIR=out
PAIRS=""
PEAK_TYPE=""

while :; do
    case $1 in

        -n|--narrow-peaks) PEAK_TYPE="-n"
        ;;

        --proper-pair-off) PROPER_PAIR="No"
        ;;

        -s|--single-end-on) SE="Yes"
        ;;

        -o|--outdir)
            ODIR=$2
            shift
        ;;

        --pairing-file)
            PAIRS=$2
            shift
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
echo SE=$SE
BAMS=$*
echo SDIR=$SDIR
echo BAMS=$BAMS
GENOME=$($SDIR/getGenomeBuildBAM.sh $1)

if [[ "$GENOME" =~ unknown ]]; then
    echo
    echo "  FATAL ERROR"
    echo "  UNKNOWN GENOME" $GENOME
    echo
    exit -1
fi

echo GENOME=$GENOME

mkdir -p $ODIR

RUNTIME="-W 59"
RUNTIMELONG="-W 359"

#if [ "" ]; then
if [ $SE = "No" ]; then

    if [ $PROPER_PAIR = "Yes" ]; then

        echo $BAMS \
            | xargs -n 1 bsub $RUNTIMELONG -o LSF.POST/ -J ${TAG}_01_POST2_$$ -R "rusage[mem=24]" \
                $SDIR/postMapBamProcessing_ChIPSeq.sh $ODIR

    else

        echo $BAMS \
            | xargs -n 1 bsub $RUNTIMELONG -o LSF.POST/ -J ${TAG}_01_POST2_$$ -R "rusage[mem=24]" \
                $SDIR/postMapBamProcessing_ChIPSeq_NoPP.sh $ODIR


    fi

else

    echo $BAMS \
        | xargs -n 1 bsub $RUNTIMELONG -o LSF.POST/ -J ${TAG}_01_POST2_$$ -R "rusage[mem=24]" \
            $SDIR/postMapBamProcessing_ChIPSeq_SE.sh $ODIR

fi

bSync ${TAG}_01_POST2_$$

ls $ODIR/*.bed.gz \
    | xargs -n 1 bsub $RUNTIMELONG -o LSF.BW/ -J ${TAG}_02_BW2_$$ -R "rusage[mem=24]" \
        $SDIR/makeBigWigFromBEDZ.sh $GENOME

bSync ${TAG}_02_BW2_$$

medianFragmentLength=$(Rscript --no-save $SDIR/getMedianFragmentLengthFromPredictDFile.R $ODIR/profiles/*.log)

echo "medianFragmentLength =" $medianFragmentLength

if [ "$PAIRS" == "" ]; then
    echo
    echo "Unpaired peak calling not yet implemented"
    echo "To run unpaired samples put 'na' in col 1"
    echo "of pairing file"
    echo
    exit
fi

Rscript --no-save $SDIR/generateMACSArgs.R $PAIRS $ODIR/*.bed.gz \
    | xargs -n 2 bsub $RUNTIMELONG -o LSF.CALLP/ -J ${TAG}_03_CALLP2_$$ -n 3 -R "rusage[mem=10]" \
        $SDIR/callPeaks_ChIPseq.sh $PEAK_TYPE $GENOME $medianFragmentLength

bSync ${TAG}_03_CALLP2_$$

#fi # DEBUG

bsub $RUNTIME -o LSF.POST/ -J ${TAG}_MergePeaks_$$ -n 3 -R "rusage[mem=10]" \
    $SDIR/mergePeaksToSAF.sh $ODIR/macs \>$ODIR/macs/macsPeaksMerged.saf

PBAMS=$(ls $ODIR/*_postProcess.bam)
bsub $RUNTIMELONG -o LSF.POST/ -J ${TAG}_Count_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_MergePeaks_$$)" \
    $SDIR/featureCounts -O -Q 10 -p -T 10 \
        -F SAF -a $ODIR/macs/macsPeaksMerged.saf \
        -o $ODIR/macs/peaks_raw_fcCounts.txt \
        $PBAMS

bsub $RUNTIME -o LSF.DESEQ/ -J ${TAG}_DESEQ_$$ -R "rusage[mem=24]" -w "post_done(${TAG}_Count_$$)" \
    Rscript --no-save $SDIR/getDESeqScaleFactors.R $ODIR/macs/peaks_raw_fcCounts.txt

bSync ${TAG}_DESEQ_$$
module unload bedtools

date >00_PIPE_DONE

cat << EOF | tee -a 00_PIPE_DONE

pipe.sh done

  Run
      ChIP-Seq/postPipe.sh out

  after checking output and LSF logs

  $ODIR=out

EOF
