#!/bin/bash
#
# genomeInfo.sh TAG FIELD
# genomeInfo.sh --tags
#
# Print one value from lib/genomes/GENOMES.tsv, the table that holds
# every genome-dependent value of the pipeline. A supported genome is a
# row of that table; nothing else in the pipeline branches on a genome
# name. Exits 1, naming the supported tags or fields, for an unknown TAG
# or FIELD.
#
# FIELD is a column name of the table (SPECIES, MACS_GSIZE, HOMER, TXDB,
# ORGDB, ADD_CHR) or GENOME_FILE, the path of lib/genomes/TAG.genome
# (autosomes, X and Y: the POST region set and the chromosome sizes).
#
# --tags prints the supported tags, space separated.
#

SDIR="$( cd "$( dirname "$0" )" && pwd )"
TABLE=$SDIR/../lib/genomes/GENOMES.tsv

if [ ! -s "$TABLE" ]; then
    echo "FATAL ERROR: genome table missing [$TABLE]" >&2
    exit 1
fi

TAGS=$(awk -F'\t' 'NR > 1 {printf "%s%s", sep, $1; sep = " "}' "$TABLE")

if [ "$#" == "1" ] && [ "$1" == "--tags" ]; then
    echo $TAGS
    exit 0
fi

if [ "$#" != "2" ]; then
    echo "usage: genomeInfo.sh TAG FIELD | genomeInfo.sh --tags" >&2
    exit 1
fi

TAG=$1
FIELD=$2

if ! awk -F'\t' -v t="$TAG" 'NR > 1 && $1 == t {found = 1} END {exit !found}' "$TABLE"; then
    echo "FATAL ERROR: unsupported genome [$TAG]; supported: $TAGS" >&2
    exit 1
fi

if [ "$FIELD" == "GENOME_FILE" ]; then
    GFILE=$(cd "$SDIR/../lib/genomes" && pwd)/$TAG.genome
    if [ ! -s "$GFILE" ]; then
        echo "FATAL ERROR: genome file missing [$GFILE]" >&2
        exit 1
    fi
    echo $GFILE
    exit 0
fi

VALUE=$(awk -F'\t' -v t="$TAG" -v f="$FIELD" '
    NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
    $1 == t && (f in col) {print $col[f]; found = 1}
    END {exit !found}' "$TABLE")

if [ "$?" != "0" ] || [ -z "$VALUE" ]; then
    FIELDS=$(head -1 "$TABLE" | cut -f2- | tr '\t' ' ')
    echo "FATAL ERROR: unknown field [$FIELD]; fields: $FIELDS GENOME_FILE" >&2
    exit 1
fi

echo $VALUE
