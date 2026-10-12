#!/bin/bash
#
# callPeaks_ChIPseq.sh [-n] GENOME FRAGLEN TBED CBED
#
# CALLP stage: macs2 callpeak on one target/control pair of clean.bed.gz
# files (a line of out/macsPairs.txt). CBED is _NoCTRL for a target with
# no control ('na' in the pairing file). Broad peaks unless -n.
#
# Output: out/macs/PREFIX/PREFIX_peaks.*, where PREFIX is
# <target>_<control>, the two sample names, or <target>____NoCTRL.
#

set -e -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

. $SDIR/venv/bin/activate
MACS=$SDIR/venv/bin/macs2

usage() {
    echo
    echo "    usage: callPeaks_ChIPseq.sh [-n] GENOME FRAGLEN TBED CBED|_NoCTRL"
    echo
    exit 1
}

CALL_BROAD_PEAKS="--broad"
while :; do
    case ${1:-} in
        -n|--narrow-peaks) CALL_BROAD_PEAKS=""
        echo
        echo "Calling narrow peaks set"
        echo
        ;;

        -*)
        echo
        echo "Invalid option [$1]"
        usage
        ;;

        *) break
        ;;

    esac
    shift
done

if [ "$#" != "4" ]; then
    usage
fi

GBUILD=$1
FRAGSIZE=$2
TBED=$3
CBED=$4

MACS_GENOME=$($SDIR/bin/genomeInfo.sh $GBUILD MACS_GSIZE) \
    || { echo "FATAL ERROR: no MACS genome size for [$GBUILD]"; exit 1; }

echo "MACS GENOME =" $MACS_GENOME

if [[ ! "$FRAGSIZE" =~ ^[1-9][0-9]*$ ]]; then
    echo "FATAL ERROR: fragment length must be an integer > 0 [$FRAGSIZE]"
    exit 1
fi

SUFFIX=_postProcess.clean.bed.gz

sampleName() {
    local name
    name=$(basename $1)
    if [[ "$name" != *$SUFFIX ]]; then
        echo "FATAL ERROR: not a POST clean BED [$1]" >&2
        return 1
    fi
    echo ${name%$SUFFIX}
}

TNAME=$(sampleName $TBED)

if [ "$CBED" == "_NoCTRL" ]; then
    CNAME="___NoCTRL"
    MACS_C_ARG=""
else
    CNAME=$(sampleName $CBED)
    MACS_C_ARG="-c $CBED"
fi

PREFIX=${TNAME}_${CNAME}
echo PREFIX=$PREFIX

ODIR=$(dirname $TBED)/macs/$PREFIX
mkdir -p $ODIR

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

#
# MACS2 writes its temp files to TMPDIR, the per-job directory on
# /localscratch that bin/runStage.sh creates.
#

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
    --outdir $ODIR \
    || { echo "FATAL ERROR: macs2 callpeak failed for $PREFIX"; exit 1; }
