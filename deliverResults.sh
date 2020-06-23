#!/bin/bash

RESDIR=$1

mkdir -p $RESDIR/chipSeq/macs
mkdir -p $RESDIR/chipSeq/bw

if [ -e "./annote" ]; then
    echo "Annotation exists"
    mkdir -p $RESDIR/chipSeq/annote
fi

rsync -rvP out/profiles/*.bw $RESDIR/chipSeq/bw
rsync -rvP out/macs $RESDIR/chipSeq
