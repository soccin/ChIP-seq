#!/bin/bash

export PYTHONPATH=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/lib/python2.7/site-packages/:$PYTHONPATH
MACS=/opt/common/CentOS_6/MACS2/MACS2-2.1.1/bin/macs2

SDIR="$( cd "$( dirname "$0" )" && pwd )"

GBUILD=$1

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

TBED=$2
if [ "$#" == "3" ]; then
    CBED=$3
else
    CBED="_NoCTRL"
fi

ODIR=$(dirname $TBED)/macs
mkdir -p $ODIR

# TDIR=/scratch/socci/_scratch_ATACSeq/$(uuidgen -t)
# mkdir -p $TDIR
# echo $TDIR, $ODIR
# export TMPDIR=$TDIR

#
# From:
#   https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/macs2/src/macs2.py
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


smooth_window=150

# ./modules/callpeak_macs2_atac.bds:    shiftsize := round( smooth_window.parseReal()/2.0 )

shiftsize=$((smooth_window / 2))

# From paper (Philip, et al, )

pval_thres=0.01

$MACS callpeak \
    -t $TDIR/cleanBED.bed.gz \
    -f BED \
    -n $PREFIX \
    -g $genome \
    -p $pval_thres \
    --nomodel \
    --shift $shiftsize \
    --extsize $smooth_window \
    --call-summits \
    --outdir $ODIR

MACS_ERROR=$?

exit $MACS_ERROR
