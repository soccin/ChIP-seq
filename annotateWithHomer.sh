#!/bin/bash
#
# annotateWithHomer.sh MACS_PEAK_FILE
#
# Annotate one MACS peak file (_peaks.broadPeak or _peaks.narrowPeak)
# with HOMER annotatePeaks.pl. Run from the analysis directory of a
# completed pipe.sh run. Writes annote/<peak file>_HOMERAnnote.txt and
# .xlsx.
#
# The genome is GENOME= in 00.RUNSTATUS.txt (detected by pipe.sh from the
# BAMs); the HOMER genome and whether to add "chr" to the contig names
# come from lib/genomes/GENOMES.tsv (b37: HOMER hg19 with "chr" added;
# hg38 and mm10: unchanged). Every peak must be on a contig of
# lib/genomes/<GENOME>.genome, so a peak file of another run or build is
# an error rather than a file of unannotated peaks.
#
# On IRIS:
#
#   find out/macs | egrep "_peaks.(broad|narrow)Peak$" \
#       | xargs -n 1 bsub -o SLURM.HOMER/ -J HOMER -W 59 -n 2 -M 16G \
#           ./ChIP-seq/annotateWithHomer.sh
#

set -eu -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

HOMERDIR=${CHIP_HOMER_DIR:-/usersoftware/core001/common/RHEL_8/homer}
HOMERBIN=$HOMERDIR/bin

if [ "$#" != "1" ]; then
    echo
    echo "  usage: annotateWithHomer.sh MACS_PEAK_FILE"
    echo
    echo "     _peaks.broadPeak or _peaks.narrowPeak; run from the analysis"
    echo "     directory. The genome is read from 00.RUNSTATUS.txt."
    echo
    exit 1
fi

PEAKFILE=$1

if [ ! -x "$HOMERBIN/annotatePeaks.pl" ]; then
    echo
    echo "    FATAL ERROR: no HOMER annotatePeaks.pl in [$HOMERBIN]"
    echo "    Install HOMER (http://homer.ucsd.edu/homer/configureHomer.pl)"
    echo "    or set CHIP_HOMER_DIR"
    echo
    exit 1
fi

if [ ! -r "$PEAKFILE" ]; then
    echo "FATAL ERROR: cannot read peak file [$PEAKFILE]"
    exit 1
fi

if [ ! -s 00.RUNSTATUS.txt ]; then
    echo "FATAL ERROR: no 00.RUNSTATUS.txt in $PWD; run from the analysis directory of a pipe.sh run"
    exit 1
fi

GENOME=$(sed -n 's/^GENOME=//p' 00.RUNSTATUS.txt | tail -1)
if [ -z "$GENOME" ]; then
    echo "FATAL ERROR: no GENOME in 00.RUNSTATUS.txt"
    exit 1
fi

HOMER_GENOME=$($SDIR/bin/genomeInfo.sh $GENOME HOMER)
ADD_CHR=$($SDIR/bin/genomeInfo.sh $GENOME ADD_CHR)
GFILE=$($SDIR/bin/genomeInfo.sh $GENOME GENOME_FILE)

if [ ! -d "$HOMERDIR/data/genomes/$HOMER_GENOME" ]; then
    echo "FATAL ERROR: HOMER genome [$HOMER_GENOME] for [$GENOME] is not installed in $HOMERDIR/data/genomes"
    exit 1
fi

echo GENOME=$GENOME HOMER_GENOME=$HOMER_GENOME ADD_CHR=$ADD_CHR
echo PEAKFILE=$PEAKFILE

#
# Every peak must lie on a contig of the run's genome file and end within
# it. Contig names alone do not tell hg38 from mm10 (both have chr1 to
# chr19); the lengths do.
#
BAD=$(awk -F'\t' 'NR == FNR {len[$1] = $2; next}
                   (!($1 in len) || $3 > len[$1]) && !($1 in seen) && n < 5 {
                       seen[$1] = 1; n++; printf "%s:%s ", $1, $3}' \
          $GFILE $PEAKFILE)
if [ -n "$BAD" ]; then
    echo "FATAL ERROR: $PEAKFILE has peaks outside the contigs of $GFILE ($GENOME): $BAD"
    exit 1
fi

export PATH=$HOMERBIN:$PATH
mkdir -p annote

TMPFILE=$(mktemp -p "$TMPDIR" homerPeaks.XXXXXX.bed)
trap 'rm -f "$TMPFILE"' EXIT

ANNOTEFILE=annote/$(basename ${PEAKFILE%.txt})_HOMERAnnote.txt

if [ "$ADD_CHR" == "yes" ]; then
    awk -F'\t' 'BEGIN {OFS = "\t"} {$1 = "chr"$1; print}' $PEAKFILE >$TMPFILE
else
    cp $PEAKFILE $TMPFILE
fi

$HOMERBIN/annotatePeaks.pl $TMPFILE $HOMER_GENOME >$ANNOTEFILE \
    || { echo "FATAL ERROR: annotatePeaks.pl failed on $PEAKFILE"; exit 1; }

Rscript --no-init-file $SDIR/postProcessHomer.R $ANNOTEFILE \
    || { echo "FATAL ERROR: postProcessHomer.R failed on $ANNOTEFILE"; exit 1; }

echo "annotateWithHomer: $ANNOTEFILE"
