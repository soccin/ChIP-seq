#!/bin/bash

export PYTHONPATH=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/lib/python2.7/site-packages/:$PYTHONPATH
MACS=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/bin/macs2

SDIR="$( cd "$( dirname "$0" )" && pwd )"

CALL_BROAD_PEAKS="--broad"
while :; do
    case $1 in
        -n|--narrow-peaks) CALL_BROAD_PEAKS=" "
        echo
        echo "Calling narrow peaks set"
        echo
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

GBUILD=$1
FRAGSIZE=$2

case $GBUILD in
    mm10)
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

echo "MACS GENEOME =" $MACS_GENOME

TBED=$3
if [ "$#" == "4" ]; then
    CBED=$4
    MACS_C_ARG="-c $CBED"
else
    CBED="___NoCTRL"
    MACS_C_ARG=" "
fi

PREFIX=$(basename ${TBED/_postProcess.*/})_$(basename ${CBED/_postProcess.*/} | sed 's/.*_s_/___/')
echo $PREFIX

ODIR=$(dirname $TBED)/macs/$PREFIX
mkdir -p $ODIR

# TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
# mkdir -p $TDIR
# echo $TDIR, $ODIR
# export TMPDIR=$TDIR

#
# From:
#   https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/macs2/src/macs2.py
#
#
# -B --SPMR ask MACS2 to generate pileup signal file of 'fragment pileup per million reads' in bedGraph format.
#

    # # ===========================================
    # # Generate narrow peaks and preliminary signal tracks
    # # ============================================

    # command = 'macs2 callpeak ' + \
    #           '-t %s -c %s ' % (experiment.name, control.name) + \
    #           '-f BED -n %s/%s ' % (peaks_dirname, prefix) + \
    #           '-g %s -p 1e-2 --nomodel --shift 0 --extsize %s --keep-dup all -B --SPMR' % (genomesize, fraglen)


    # # ===========================================
    # # Generate Broad and Gapped Peaks
    # # ============================================

    # command = 'macs2 callpeak ' + \
    #           '-t %s -c %s ' % (experiment.name, control.name) + \
    #           '-f BED -n %s/%s ' % (peaks_dirname, prefix) + \
    #           '-g %s -p 1e-2 --broad --nomodel --shift 0 --extsize %s --keep-dup all' % (genomesize, fraglen)


pval_thres=0.01

$MACS callpeak \
    -t $TBED \
    $MACS_C_ARG \
    -f BED \
    -n $PREFIX \
    -g $MACS_GENOME \
    -p $pval_thres \
    --nomodel \
    --shift 0 \
    --extsize $FRAGSIZE \
    $CALL_BROAD_PEAKS \
    --outdir $ODIR

MACS_ERROR=$?

exit $MACS_ERROR
