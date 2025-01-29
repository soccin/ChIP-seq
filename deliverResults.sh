#!/bin/bash

RESDIR=$1

if [ "$RESDIR" == "" ]; then
    echo
    echo "  usage: deliverResults.sh RESDIR"
    echo "    RESDIR = .../pi/invest/projNo/r_###"
    echo
    exit
fi

set -eu

echo ""
echo "Need sudo to chmod output folder"
echo ""

sudo chmod g+ws $RESDIR
mkdir -p $RESDIR/chipSeq/macs
mkdir -p $RESDIR/chipSeq/bw
mkdir -p $RESDIR/chipSeq/qc
mkdir -p $RESDIR/chipSeq/pipeline_info

if [ -e "out/annote" ]; then
    echo "Annotation exists"
    mkdir -p $RESDIR/chipSeq/annote
    rsync -rvP out/annote $RESDIR/chipSeq
fi

rsync -rvP out/profiles/*.bw $RESDIR/chipSeq/bw
rsync -rvP out/macs $RESDIR/chipSeq
rsync -rvP out/qc $RESDIR/chipSeq

if [ -e "out/diff" ]; then
    echo "Diff Analysis Exists"
    mkdir -p $RESDIR/chipSeq/diff
    rsync -rvP out/diff $RESDIR/chipSeq
fi

cp out/pipeline_info/* $RESDIR/chipSeq/pipeline_info

PROJNO=$(echo $RESDIR | tr '/' '\n' | fgrep Proj_ | sed 's/Proj_//')

cat <<EOF

ChIPSeq project $PROJNO results ready

The outputs for project Proj_$PROJNO are ready.

You can access them on the BIC server at:

	https://bicdelivery.mskcc.org/project/$PROJNO/chipseq/r_001

<<COPY_PASTE ChIPSeq/docs/RESULTS.md>>

EOF
