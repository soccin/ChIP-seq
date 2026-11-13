#!/bin/bash
#
# 00.SETUP.sh
#
# Builds venv/ for MACS2 in the directory that holds this script (the
# ChIP-seq clone), whatever the current directory. Run it once; do not
# source it. From the analysis directory:
#
#   ./ChIP-seq/00.SETUP.sh
#
# MACS2 2.2.9.1 builds under the python 3.10 on the IRIS PATH with numpy
# pinned to 1.23.0 (the recipe ATAC-seq uses on IRIS).
#
# Stops if venv/ already exists. If a step fails, venv/ is left
# incomplete; remove it and run again.
#

# Sourced, the exit calls below would close the calling shell.
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
    echo "FATAL ERROR: run 00.SETUP.sh, do not source it" >&2
    return 1
fi

set -eu -o pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"
cd "$SDIR"

if [ -e venv ]; then
    echo "FATAL ERROR: $SDIR/venv already exists; remove it to rebuild" >&2
    exit 1
fi

INCOMPLETE="$SDIR/venv is incomplete; remove it and run again"

python3 -m venv venv \
    || { echo "FATAL ERROR: python3 -m venv failed; $INCOMPLETE" >&2; exit 1; }
venv/bin/python -m pip install --upgrade pip \
    || { echo "FATAL ERROR: pip upgrade failed; $INCOMPLETE" >&2; exit 1; }
venv/bin/python -m pip install numpy==1.23.0 \
    || { echo "FATAL ERROR: numpy 1.23.0 install failed; $INCOMPLETE" >&2; exit 1; }
venv/bin/python -m pip install MACS2==2.2.9.1 \
    || { echo "FATAL ERROR: MACS2 2.2.9.1 install failed; $INCOMPLETE" >&2; exit 1; }
