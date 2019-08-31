#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

ODIR=$1
IBAM=$2
OBAM=$ODIR/$(basename ${IBAM/.bam/_postProcess.bam})

#TDIR=/scratch/socci/_scratch_ChIPSeq/$(uuidgen -t)
TDIR=_scratch_ChIPSeq/$(uuidgen -t)
mkdir -p $TDIR
echo TDIR=$TDIR

echo ODIR=$ODIR
echo OBAM=$OBAM


#
# Generate a BED file with the chromosomes to keep
# i.e, remove non standard ones and MT
#
samtools view -H $IBAM \
    | fgrep SQ \
    | egrep "SN:([1-9XY]|chr[1-9XY])" \
    | egrep -v "_" \
    | cut -f2-3 \
    | sed 's/..://g' \
    | awk '{print $1"\t0\t"$2}' \
    > $TDIR/regionsToKeep_$$


#
# From ENCODE
#    https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/filter_qc/src/filter_qc.py
#

MAPQ=30

# f 2 ==> proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
## F 3852 (1804+2048) 2048==Supplementary Alignment
samtools view -q $MAPQ -F 3852 -L $TDIR/regionsToKeep_$$ $IBAM -u >$TDIR/step1.bam

picardV2 SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true

#
# Create BED version but keep filter BAM also
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | gzip -nc >${OBAM/.bam/.clean.bed.gz}

#rm -rf $TDIR
