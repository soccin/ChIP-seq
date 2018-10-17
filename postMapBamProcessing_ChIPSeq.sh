#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

IBAM=$1

if [ "$#" == "1" ]; then
    OBAM=${IBAM/.bam/_postProcess.bam}
    OBAM=$(basename $OBAM)
    echo $OBAM
else
    OBAM=$2
fi

TDIR=/scratch/socci/_scratch_ChIPSeq/$(uuidgen -t)
#TDIR=_scratch_ChIPSeq/$(uuidgen -t)
mkdir -p $TDIR
echo $TDIR

#
# Generate a BED file with the chromosomes to keep
# i.e, remove non standard ones and MT
#
samtools view -H $IBAM \
    | fgrep SQ \
    | egrep "SN:([1-9XY]|chr[1-9XY])" \
    | cut -f2-3 \
    | sed 's/..://g' \
    | awk '{print $1"\t0\t"$2}' \
    > $TDIR/regionsToKeep_$$


# f 3 ==> paired, proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
samtools view -q 10 -f 3 -F 1804 -L $TDIR/regionsToKeep_$$ $IBAM -u >$TDIR/step1.bam
picardV2 SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true

#
# Create BED version but keep filter BAM also
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | gzip -nc >${OBAM/.bam/.clean.bed.gz}

#rm -rf $TDIR
