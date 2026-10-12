#!/bin/bash
#
# postMapBamProcessing_ChIPSeq_SE.sh GENOME SAMPLE BAM
#
# POST stage for pipe.sh -s (single-end BAMs): the ENCODE3 filters that
# apply to single reads (MAPQ, flags, contigs), then a coordinate sort.
# Writes, under out/:
#
#   SAMPLE_postProcess.bam, .bai and .bam.bai (relative link to the .bai)
#   SAMPLE_postProcess.clean.bed.gz
#

set -e -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"

. $SDIR/bin/loadTools.sh

if [ "$#" != "3" ]; then
    echo
    echo "    usage: postMapBamProcessing_ChIPSeq_SE.sh GENOME SAMPLE BAM"
    echo
    exit 1
fi

GENOME=$1
SAMPLE=$2
IBAM=$3

ODIR=out
OBAM=$ODIR/${SAMPLE}_postProcess.bam

GFILE=$($SDIR/bin/genomeInfo.sh $GENOME GENOME_FILE) \
    || { echo "FATAL ERROR: no genome file for [$GENOME]"; exit 1; }

TDIR=${CHIP_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ChIPSeq}/$(uuidgen -t)
mkdir -p $TDIR
trap 'rm -rf "$TDIR"' EXIT

echo GENOME=$GENOME SAMPLE=$SAMPLE
echo IBAM=$IBAM
echo OBAM=$OBAM
echo TDIR=$TDIR

#
# Regions to keep: the contigs of the genome file, whole. Removes MT and
# every unplaced, random, decoy and viral contig.
#
awk 'BEGIN {OFS = "\t"} {print $1, 0, $2}' $GFILE >$TDIR/regionsToKeep.bed \
    || { echo "FATAL ERROR: cannot write regions from $GFILE"; exit 1; }

#
# From ENCODE
#    https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/filter_qc/src/filter_qc.py
#

MAPQ=30

# F 1804 ==> unmapped, mate unmapped, not primary, fails QC, duplicate
## F 3852 (1804+2048) 2048==Supplementary Alignment
samtools view -q $MAPQ -F 3852 -L $TDIR/regionsToKeep.bed $IBAM -u >$TDIR/step1.bam \
    || { echo "FATAL ERROR: samtools filter failed on $IBAM"; exit 1; }

$SDIR/bin/picardV2 SortSam I=$TDIR/step1.bam O=$OBAM SO=coordinate MAX_RECORDS_IN_RAM=5000000 CREATE_INDEX=true \
    || { echo "FATAL ERROR: picard SortSam coordinate failed"; exit 1; }
rm $TDIR/step1.bam

if [ ! -s "${OBAM%.bam}.bai" ]; then
    echo "FATAL ERROR: no index written for $OBAM"
    exit 1
fi
ln -sfn $(basename ${OBAM%.bam}.bai) $OBAM.bai

#
# Create BED version but keep filter BAM also
#

samtools view -b $OBAM \
    | bedtools bamtobed -i - \
    | gzip -nc >${OBAM%.bam}.clean.bed.gz \
    || { echo "FATAL ERROR: BED conversion failed for $OBAM"; exit 1; }
