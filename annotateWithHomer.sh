#/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

HOMERBIN=/ifs/work/socci/tools/HOMER/bin
export PATH=$HOMERBIN:$PATH
mkdir -p annote

macsNarrowPeakFile=$1

wdir=$(dirname $macsNarrowPeakFile)
TMPFILE=$wdir/tmpPeaks.bed_$$
echo TMPFILE=$TMPFILE

ANNOTEFILE=annote/$(basename ${macsNarrowPeakFile/.txt/}_HOMERAnnote.txt)

cat $macsNarrowPeakFile | awk '{print "chr"$0}' >$TMPFILE
$HOMERBIN/annotatePeaks.pl $TMPFILE hg19 >$ANNOTEFILE
echo rm $TMPFILE

Rscript --no-save $SDIR/postProcessHomer.R $ANNOTEFILE
