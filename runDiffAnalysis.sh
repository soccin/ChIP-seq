#!/bin/bash
SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ ! -e "annote/macsPeaksMerged.bed_HOMERAnnote.txt" ]; then

	cat out/macs/macsPeaksMerged.saf \
		| awk '{print $2,$3,$4,$1}' \
	    | tr ' ' '\t' >macsPeaksMerged.bed

	echo
	echo Only works for human
	echo

	./ChIP-seq/annotateWithHomer.sh macsPeaksMerged.bed
fi

Rscript ChIP-seq/diffAnalysisChIP.R

mkdir -pv out/diff
mv -v *__diffAnalysis__* out/diff


