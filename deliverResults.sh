#!/bin/bash

RESDIR=$1

if [ "$RESDIR" == "" ]; then
    echo
    echo "  usage: deliverResults.sh RESDIR"
    echo "    RESDIR = .../pi/invest/projNo/r_###"
    echo
    exit
fi

mkdir -p $RESDIR/chipSeq/macs
mkdir -p $RESDIR/chipSeq/bw

if [ -e "./annote" ]; then
    echo "Annotation exists"
    mkdir -p $RESDIR/chipSeq/annote
    rsync -rvP annote $RESDIR/chipSeq
fi

rsync -rvP out/profiles/*.bw $RESDIR/chipSeq/bw
rsync -rvP out/macs $RESDIR/chipSeq
