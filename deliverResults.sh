#!/bin/bash

RESDIR=$1

if [ "$RESDIR" == "" ]; then
    echo
    echo "  usage: deliverResults.sh RESDIR"
    echo "    RESDIR = .../pi/invest/projNo/r_###"
    echo
    exit
fi

echo ""
echo "Need sudo to chmod output folder"
echo ""

sudo chmod g+ws $RESDIR
mkdir -p $RESDIR/chipSeq/macs
mkdir -p $RESDIR/chipSeq/bw
mkdir -p $RESDIR/chipSeq/qc

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

