#/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

HOMERBIN=/lila/data/bicgrp/socci/work/homer/bin

if [ ! -e "$HOMERBIN" ]; then
    echo -e "\n\nNeed to install HOMER\n"
    echo -e "Start with  wget http://homer.ucsd.edu/homer/configureHomer.pl"
    echo -e "and set HOMERBIN [current:${HOMERBIN}]"
    echo
    exit
fi

export PATH=$HOMERBIN:$PATH
mkdir -p annote

macsNarrowPeakFile=$1

if [ "$macsNarrowPeakFile" == "" ]; then
    echo
    echo "  usage: annotateWithHomer.sh MACS_PEAK_FILE"
    echo
    echo "     Can use either _peaks.broadPeak or _peaks.narrowPeak "
    echo "     N.B.: Hardcoded to HG19"
    echo
    exit
fi

echo
echo
echo "##########################################"
echo This is hardcoded to HG19
echo
echo

wdir=$(dirname $macsNarrowPeakFile)
TMPFILE=$wdir/tmpPeaks.bed_$$
echo TMPFILE=$TMPFILE

ANNOTEFILE=annote/$(basename ${macsNarrowPeakFile/.txt/}_HOMERAnnote.txt)

cat $macsNarrowPeakFile | awk '{print "chr"$0}' >$TMPFILE
$HOMERBIN/annotatePeaks.pl $TMPFILE hg19 >$ANNOTEFILE
echo rm $TMPFILE

Rscript --no-save $SDIR/postProcessHomer.R $ANNOTEFILE
