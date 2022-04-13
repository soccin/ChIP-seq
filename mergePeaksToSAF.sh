#!/bin/bash

MACS_OUTDIR=$1

find $MACS_OUTDIR \
    | egrep "narrowPeak$|broadPeak$" \
    | xargs sort -S 16g -k1,1V -k2,2n \
    | bedtools merge -i - -d 500 \
    | awk '{print "Peak_"++s,$0,"+"}' | tr ' ' '\t'

