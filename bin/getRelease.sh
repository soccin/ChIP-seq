#!/bin/bash

#
# getRelease.sh -- print the release of this checkout, vX.Y.Z, from
# the "Version: X.Y.Z" line of VERSION.md
#
# The release is the client-facing version: it is set by hand in
# VERSION.md when a release is made, never derived from git (pipe.sh
# records `git describe` as VERSION= in 00.RUNSTATUS.txt).
# deliverResults.sh puts the release in the email and the doc link, so
# the tag v<X.Y.Z> must be pushed.
#

set -euo pipefail

SDIR="$( cd "$( dirname "$0" )" && pwd )"
VFILE=$(dirname "$SDIR")/VERSION.md

RELEASE=$(sed -n 's/^Version: *//p' "$VFILE" 2>/dev/null || true)
if [[ ! "$RELEASE" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "getRelease.sh: no 'Version: X.Y.Z' line in $VFILE" >&2
    exit 1
fi

echo "v$RELEASE"
