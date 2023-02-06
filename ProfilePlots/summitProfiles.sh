#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

SUMMITS=$1

BASE=$(basename ${SUMMITS/.bed/})

TARGET=$(echo $BASE | sed 's/.*_s_/s_/' | sed 's/____.*//')
OUTDIR=$(realpath $(dirname $SUMMITS | sed 's/macs.*//'))
TARGET_BAM=$(find $OUTDIR -name '*.bam' | fgrep ${TARGET}_postProcess.bam)

Rscript $SDIR/genSummitProfileBed.R $SUMMITS

PROFILE_BED=$(basename ${SUMMITS/.bed/__Profile.bed})

bedtools multicov \
    -bed $PROFILE_BED \
    -bams $TARGET_BAM \
    > ${PROFILE_BED/.bed/}____${TARGET}.cov.txt



