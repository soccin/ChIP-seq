#!/bin/bash
#
# postPipe.sh
#
# Run from the analysis directory after pipe.sh. Refuses to run unless
# bin/checkRun.sh says the run worked. Then:
#
#   - moves out/*.bam, .bai into out/bam/ and out/*.bed.gz into out/bed/
#     (once; the count checks still find the BAMs there)
#   - QC stage 1 (qc_ChIPSeq_01.R): peak counts, significant peaks and
#     volcano plots per MACS call, into out/qc/
#   - QC stage 2 (qc_ChIPSeq_02.R) if manifest.txt (SampleID<TAB>Group, no
#     header) exists; otherwise prints the command
#
# Output is always out/ (pipe.sh -o was removed in v1.0.0), so this takes
# no arguments.
#

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long --tags --dirty)
PIPENAME="ChIP-Seq"

function usage {
    echo
    echo "usage: $PIPENAME/postPipe.sh"
    echo "version=$SCRIPT_VERSION"
    echo
    echo "    Run from the analysis directory of a completed pipe.sh run."
    echo "    Output is always out/; there are no arguments."
    echo
    exit 1
}

if [ "$#" != "0" ]; then
    usage
fi

ODIR=out

if ! $SDIR/bin/checkRun.sh; then
    echo
    echo "    FATAL ERROR: bin/checkRun.sh does not pass this run; postPipe.sh not run"
    echo
    exit 1
fi

if [ ! -s out/pipeline_info/samples.tsv ]; then
    echo "    FATAL ERROR: no out/pipeline_info/samples.tsv"
    exit 1
fi

if [ ! -e "$ODIR/bam" ]; then
    mkdir $ODIR/bam
    if compgen -G "$ODIR/*.ba?" >/dev/null; then
        mv $ODIR/*.ba? $ODIR/bam
    fi

    mkdir $ODIR/bed
    if compgen -G "$ODIR/*.bed.gz" >/dev/null; then
        mv $ODIR/*.bed.gz $ODIR/bed
    fi
fi

if [ ! -d "$ODIR/macs" ]; then
    echo
    echo "    No out/macs (run without --pairing-file); no peak QC to do."
    echo
    exit 0
fi

Rscript --no-init-file $SDIR/qc_ChIPSeq_01.R

mkdir -p $ODIR/qc/peaks
for PATTERN in "*___sigPeaks*" "*___Volcano*"; do
    if compgen -G "$PATTERN" >/dev/null; then
        mv $PATTERN $ODIR/qc/peaks
    fi
done
cp qcChIPSeq_* $ODIR/qc/

date >01_POST_DONE
cat << EOF | tee -a 01_POST_DONE

Should create a manifest/group file "manifest.txt" for
stage 2 qc or use _sample_grouping.txt

and then run
    Rscript --no-init-file $SDIR/qc_ChIPSeq_02.R manifest.txt
    mv qcChIPSeq_*.pdf qcChIPSeq_*.xlsx $ODIR/qc

EOF

if [ -e manifest.txt ]; then
    Rscript --no-init-file $SDIR/qc_ChIPSeq_02.R manifest.txt
    cp qcChIPSeq_* $ODIR/qc/

cat << EOF | tee -a 01_POST_DONE

Maybe run HOMER annotation (the genome is taken from 00.RUNSTATUS.txt)

    find out/macs | egrep "_peaks.(broad|narrow)Peak$" \\
        | xargs -n 1 bsub -o SLURM.HOMER/ -J HOMER -W 59 -n 2 -M 16G \\
            $SDIR/annotateWithHomer.sh

EOF

fi
