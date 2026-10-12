#!/bin/bash

if [ -e /etc/system-release ]; then
    OSSTR=$(cat /etc/system-release)
    OSVER=$(echo $OSSTR | perl -ne 'm/ ((\d+)\.(\d+)(|.\d+)) /; print $1')
    echo $OSVER
else
    echo "MISSING /etc/system-release"
    exit 1
fi
