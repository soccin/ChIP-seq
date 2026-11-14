#!/bin/bash
#
# postPipe.sh
#
# Run from the analysis directory after pipe.sh, and again after editing
# manifest.csv. Refuses to run unless bin/checkRun.sh says the run
# worked. Then:
#
#   - moves out/*.bam, .bai into out/bam/ and out/*.bed.gz into out/bed/
#     (once; the count checks still find the BAMs there)
#   - writes manifest.csv if there is none: header MapID,SampleID,Group,
#     one row per sample of out/pipeline_info/samples.tsv, each sample
#     its own group. Edit Group (and SampleID, the name in the plots) and
#     run postPipe.sh again. R/diffAnalysisPairwise.R reads the same file.
#   - QC stage 1 (qc_ChIPSeq_01.R): peak counts into out/qc/, significant
#     peaks and volcano plots per MACS call into out/qc/peaks/
#   - QC stage 2 (qc_ChIPSeq_02.R manifest.csv): fragments in peaks per
#     sample and a PCA of the peak counts colored by Group, into out/qc/
#   - 01_POST_DONE: what was done and the next steps (HOMER annotation,
#     differential analysis)
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
    echo "    Output is always out/; there are no arguments. Writes"
    echo "    manifest.csv (each sample its own group) if there is none;"
    echo "    edit its Group column and run again to group the PCA."
    echo
    exit 1
}

if [ "$#" != "0" ]; then
    usage
fi

ODIR=out
SAMPLES_TSV=$ODIR/pipeline_info/samples.tsv
MANIFEST=manifest.csv

if ! $SDIR/bin/checkRun.sh; then
    echo
    echo "    FATAL ERROR: bin/checkRun.sh does not pass this run; postPipe.sh not run"
    echo
    exit 1
fi

if [ ! -s $SAMPLES_TSV ]; then
    echo "    FATAL ERROR: no $SAMPLES_TSV"
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

#
# Up to v1.0.0 QC stage 2 read manifest.txt (SampleID<TAB>Group, no
# header). Its groups would be ignored, so stop until it is converted.
#
if [ -e manifest.txt ] && [ ! -e "$MANIFEST" ]; then
    echo
    echo "    FATAL ERROR: manifest.txt is no longer read. Write its groups into $MANIFEST:"
    echo "        a CSV with the header MapID,SampleID,Group and every sample of the run"
    echo "        (example: $SDIR/assets/manifest_example.csv)"
    echo
    exit 1
fi

WROTE_MANIFEST=No
if [ ! -e "$MANIFEST" ]; then
    awk -F'\t' '
        NR == 1 { if ($1 != "SAMPLE") exit 1; print "MapID,SampleID,Group"; next }
        { print $1 "," $1 "," $1 }
    ' $SAMPLES_TSV >$MANIFEST.tmp \
        || { echo "    FATAL ERROR: cannot write $MANIFEST from $SAMPLES_TSV"; exit 1; }
    mv $MANIFEST.tmp $MANIFEST
    WROTE_MANIFEST=Yes
    echo
    echo "    Wrote $MANIFEST: each sample its own group"
    echo
fi

Rscript --no-init-file $SDIR/qc_ChIPSeq_01.R \
    || { echo "    FATAL ERROR: QC stage 1 (qc_ChIPSeq_01.R) failed"; exit 1; }

Rscript --no-init-file $SDIR/qc_ChIPSeq_02.R $MANIFEST \
    || { echo "    FATAL ERROR: QC stage 2 (qc_ChIPSeq_02.R $MANIFEST) failed"; exit 1; }

{
date
cat << EOF

postPipe.sh $SCRIPT_VERSION
QC reports in $ODIR/qc/; PCA groups from $MANIFEST
EOF

if [ "$WROTE_MANIFEST" == "Yes" ]; then
cat << EOF

$MANIFEST was written by this run, with each sample its own group. To
group samples in the PCA, edit its Group column (SampleID is the name in
the plots; MapID must stay the sample name) and run postPipe.sh again.
Groups starting with EXC are left out of the PCA and the differential
analysis.
EOF
fi

cat << EOF

Maybe run HOMER annotation (the genome is taken from 00.RUNSTATUS.txt)

    find out/macs | egrep "_peaks.(broad|narrow)Peak$" \\
        | xargs -n 1 bsub -o SLURM.HOMER/ -J HOMER -W 59 -n 2 -M 16G \\
            $SDIR/annotateWithHomer.sh

Differential analysis, with the groups of $MANIFEST and a comparisons
file (example: $SDIR/assets/comps_example.csv)

    Rscript --no-init-file $SDIR/R/diffAnalysisPairwise.R $MANIFEST comps.csv

EOF
} | tee 01_POST_DONE
