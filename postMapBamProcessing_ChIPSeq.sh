#!/bin/bash

set -e 

module load samtools

SDIR="$( cd "$( dirname "$0" )" && pwd )"

ODIR=$1
IBAM=$2
OBAM=$ODIR/$(basename ${IBAM/.bam/_postProcess.bam})

#TDIR=/scratch/socci/_scratch_ChIPSeq/$(uuidgen -t)
TDIR=_scratch_ChIPSeq/$(uuidgen -t)
mkdir -p $TDIR
echo "[12]" $TDIR
echo IBAM=$IBAM >> $TDIR/postMapLog
echo IBAM=$IBAM

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
# REMOVE this we do not want to force proper pairs for this one f 2 ==> proper pair
# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
## F 3852 (1804+2048) 2048==Supplementary Alignment
samtools view -q $MAPQ -f 2 -F 3852 -L $TDIR/regionsToKeep_$$ $IBAM -u >$TDIR/step1.bam

# From ENCODE
# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
#
# fill in mate coordinates, ISIZE and mate-related flags
# fixmate requires name-sorted alignment; -r removes secondary and
# unmapped (redundant here because already done above?)

picardV2 SortSam I=$TDIR/step1.bam O=$TDIR/step2.bam SO=queryname MAX_RECORDS_IN_RAM=5000000
samtools fixmate -r $TDIR/step2.bam - >$TDIR/step2b.bam
samtools view -F 1804 -f 2 -u $TDIR/step2b.bam >$TDIR/step3.bam
picardV2 SortSam I=$TDIR/step3.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true
ln -s $(realpath ${OBAM/.bam/.bai}) ${OBAM/.bam/.bam.bai}

#
# Create BED version but keep filter BAM also
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | gzip -nc >${OBAM/.bam/.clean.bed.gz}

#rm -rf $TDIR
