#!/bin/bash

#
# deliverResults.sh -- deliver a finished run to the BIC delivery area
#
# usage:
#   deliverResults.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]
#   deliverResults.sh [-n] -d|--default
#
# Run from the analysis directory, after the run has COMPLETED
# (bin/checkRun.sh exits 0) and postPipe.sh has run. It
#
#   - writes project.yaml for the bicdelivery permissions database,
#     from the run and the first README.txt or REQUEST.txt found going
#     up from here to Proj_<no> (bin/readme2yaml.R);
#   - copies to DELIVERY/r_NNN/chipSeq/
#         out/profiles/*.bw     -> bw/
#         out/macs, out/qc      -> macs/, qc/  (runs with a pairing file)
#         out/pipeline_info     -> pipeline_info/
#         out/annote            -> annote/     (annotateWithHomer.sh), if any
#         out/diff              -> diff/       (R/diffAnalysisPairwise.R), if any
#         docs/output.pdf, .html of this checkout -> docs/
#     out/bam, out/bed and the other files of out/ are not delivered;
#   - writes the delivery email to deliveryEmail_<no>_r_NNN_<yymmdd>.txt,
#     with the release of VERSION.md of this checkout;
#   - registers the run with the permissions database.
#
#   -d, --default   deliver to DROOT/<pi>/<invest>/Proj_<no>, taken
#                   from the path .../Users/<Pi>/<Invest>/Proj_<no>/...
#                   of the analysis directory (lowercased but for the
#                   project)
#   -n, --dry-run   print the commands that write to the delivery area
#                   or the database instead of running them. The email
#                   and project.yaml are still written here.
#
# Without a trailing r_NNN the run folder is the one after the highest
# r_NNN in the delivery folder, else r_001. An existing r_NNN can be
# given to add files to it, e.g. a differential analysis made later;
# nothing in it is deleted.
#
# DROOT is /data1/core002/res/bic/results. Every write to the delivery
# area runs as the service account svc_core001_bot01 (dzdo).
#
# A port of ATAC-seq's deliverResults.sh (v2.1.0), with the rsync
# options of adagio's bin/deliver.sh.
#

set -euo pipefail

#
# The physical path: the docs are copied from here by the service
# account, which cannot enter a home directory that a symlink in the
# path given may pass through.
#
SDIR="$( cd "$( dirname "$0" )" && pwd -P )"

BOT_USER=svc_core001_bot01
DROOT=/data1/core002/res/bic/results
BIC_DELIVERY=$HOME/Code/BIC/Delivery/Version2j

#
# Every delivered file and folder gets this group (rsync --chown).
# rsync -a copies the source group, which the service account cannot
# set when it is not in it.
#
DGROUP=grp_hpc_core001

#
# Options of every rsync to the delivery area: the group above; every
# folder group-writable with setgid (rsync -a copies the source mode,
# and out/qc is made without them); and no macOS metadata files
# ("._<name>", written next to a file opened from a Mac).
#
RSYNC=(rsync -aP --chown=:"$DGROUP" --chmod=Dg+ws --exclude="._*")

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

WARNINGS=()
warn() {
    echo "    WARNING: $*"
    WARNINGS+=("$*")
}

#
# The physical path, as R's getwd() gives it, so that PROJ is the
# project the R reports were named by: the first Proj_* component
# (get_project_number() in tools.R, R/diffAnalysisPairwise.R).
#
CWD=$(pwd -P)

PROJ=""
IFS=/ read -r -a PARTS <<< "$CWD"
for p in "${PARTS[@]}"; do
    if [[ "$p" == Proj_* ]]; then
        PROJ=$p
        break
    fi
done

#
# Default delivery path from .../Users/Aa/BBB/Proj_nnnnn/...
# -> aaa/bbb/Proj_nnnnn; the third component must be the project.
#
DPATH=""
if [[ "$CWD" == *"/Users/"* && -n "$PROJ" ]]; then
    rest="${CWD#*/Users/}"
    IFS=/ read -r u1 u2 proj _ <<< "$rest"
    if [[ -n "$u1" && -n "$u2" && "$proj" == "$PROJ" ]]; then
        DPATH="${u1,,}/${u2,,}/$proj"
    fi
fi

usage() {
    echo
    echo "usage: deliverResults.sh [-n] /path/to/delivery/Proj_<no>[/r_NNN]"
    echo "       deliverResults.sh [-n] -d|--default"
    echo
    echo "    -n, --dry-run   print the delivery and database commands, do not run them"
    echo
    echo "    default: ${DPATH:+$DROOT/}${DPATH:-(none)}"
    echo
}

ODIR=""
USE_DEFAULT=0
DRYRUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--default)
            USE_DEFAULT=1
            ;;
        -n|--dry-run)
            DRYRUN=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
        *)
            if [[ -n "$ODIR" ]]; then
                echo "Only one delivery folder may be given"
                usage
                exit 1
            fi
            ODIR="$1"
            ;;
    esac
    shift
done

if [[ "$USE_DEFAULT" -eq 1 && -n "$ODIR" ]]; then
    echo "Give a delivery folder or -d, not both"
    usage
    exit 1
fi

if [[ "$USE_DEFAULT" -eq 1 ]]; then
    if [[ -z "$DPATH" ]]; then
        fatal "no default delivery path for $CWD (not .../Users/<Pi>/<Invest>/Proj_<no>/...)"
    fi
    ODIR="$DROOT/$DPATH"
fi

if [[ -z "$ODIR" ]]; then
    usage
    exit 1
fi

asbot() {
    if [[ "$DRYRUN" -eq 1 ]]; then
        echo "DRYRUN: dzdo -u $BOT_USER ${*@Q}"
    else
        dzdo -u "$BOT_USER" "$@"
    fi
}

##
# The run

STATUSF=00.RUNSTATUS.txt
if [ ! -e "$STATUSF" ]; then
    fatal "no $STATUSF; run deliverResults.sh from the analysis directory"
fi

if ! "$SDIR/bin/checkRun.sh"; then
    fatal "the run has not COMPLETED (bin/checkRun.sh above); nothing delivered"
fi
echo

getKey() {
    sed -n "s/^$1=//p" "$STATUSF" | tail -1
}

GENOME=$(getKey GENOME)
if [[ -z "$GENOME" ]]; then
    fatal "no GENOME= in $STATUSF"
fi

if [[ -z "$PROJ" ]]; then
    fatal "no Proj_* directory in $CWD; cannot name the project"
fi
PROJNO=${PROJ#Proj_}

#
# The email and its doc link name the release of VERSION.md of this
# checkout, whose docs are delivered; it changes only when a release
# is made, never with git.
#
if ! VERSION=$("$SDIR/bin/getRelease.sh"); then
    fatal "no release in $SDIR/VERSION.md"
fi

DOCS=("$SDIR/docs/output.pdf" "$SDIR/docs/output.html")
for f in "${DOCS[@]}"; do
    if [ ! -s "$f" ]; then
        fatal "missing $f; render it with make -C $SDIR/docs"
    fi
done

##
# What is delivered. Every run has bigWigs and pipeline_info; a run
# with a pairing file has out/macs and, once postPipe.sh has run,
# out/qc.

OUT=out

shopt -s nullglob
BW=("$OUT"/profiles/*.bw)
shopt -u nullglob
if [[ ${#BW[@]} -eq 0 ]]; then
    fatal "no bigWig files in $OUT/profiles"
fi

if [ ! -s "$OUT/pipeline_info/cmd.sh.log" ]; then
    fatal "no $OUT/pipeline_info/cmd.sh.log"
fi
FOLDERS=("$OUT/pipeline_info")

HAVE_PEAKS=0
if [ -d "$OUT/macs" ]; then
    HAVE_PEAKS=1
    QCBASE=$OUT/qc/qcChIPSeq_$PROJ
    for f in "$QCBASE.xlsx" "$QCBASE.pdf"; do
        if [ ! -s "$f" ]; then
            fatal "no $f; run postPipe.sh before delivering"
        fi
    done
    FOLDERS+=("$OUT/macs" "$OUT/qc")
else
    warn "no $OUT/macs (run without a pairing file); only the bigWigs are delivered"
fi

#
# HOMER annotation (annotateWithHomer.sh): one .txt and one .xlsx per
# MACS peak file.
#
HAVE_ANNOTE=0
if [ -d "$OUT/annote" ]; then
    shopt -s nullglob
    ANNOTE_TXT=("$OUT"/annote/*_HOMERAnnote.txt)
    ANNOTE_XLSX=("$OUT"/annote/*_HOMERAnnote.xlsx)
    PEAKFILES=("$OUT"/macs/*/*_peaks.broadPeak "$OUT"/macs/*/*_peaks.narrowPeak)
    shopt -u nullglob
    if [[ ${#ANNOTE_XLSX[@]} -eq 0 ]]; then
        fatal "$OUT/annote has no *_HOMERAnnote.xlsx; finish annotateWithHomer.sh or remove the folder"
    fi
    if [[ ${#ANNOTE_TXT[@]} -ne ${#ANNOTE_XLSX[@]} ]]; then
        warn "${#ANNOTE_TXT[@]} HOMERAnnote.txt but ${#ANNOTE_XLSX[@]} HOMERAnnote.xlsx files in $OUT/annote"
    fi
    if [[ ${#ANNOTE_XLSX[@]} -ne ${#PEAKFILES[@]} ]]; then
        warn "${#ANNOTE_XLSX[@]} HOMERAnnote.xlsx files for ${#PEAKFILES[@]} MACS peak files"
    fi
    HAVE_ANNOTE=1
    FOLDERS+=("$OUT/annote")
fi

#
# Differential analysis in out/diff: every RUNTAG of this project,
# untagged (Proj_<no>__...) or tagged (Proj_<no>_<tag>_...). They are
# delivered with the folder; the lists only check them and pick the
# email section.
#
DIFF=$OUT/diff

shopt -s nullglob
DIFF_XLSX=("$DIFF/${PROJ}"_*_DiffPeaksEdgeR_V3.xlsx)
DIFF_PDF=("$DIFF/${PROJ}"_*_DiffPeaks_V3.pdf)
shopt -u nullglob
DIFF_FILES=("${DIFF_XLSX[@]}" "${DIFF_PDF[@]}")

if [ -d "$DIFF" ]; then
    if [[ ${#DIFF_FILES[@]} -eq 0 ]]; then
        fatal "$DIFF has no ${PROJ}_*_DiffPeaks files; finish R/diffAnalysisPairwise.R or remove the folder"
    fi
    FOLDERS+=("$DIFF")
fi
if [[ ${#DIFF_XLSX[@]} -ne ${#DIFF_PDF[@]} ]]; then
    warn "${#DIFF_XLSX[@]} DiffPeaksEdgeR_V3.xlsx but ${#DIFF_PDF[@]} DiffPeaks_V3.pdf files"
fi
for f in "${DIFF_FILES[@]}"; do
    if [[ "$f" != "$DIFF/${PROJ}__DiffPeaksEdgeR_V3.xlsx" && "$f" != "$DIFF/${PROJ}__DiffPeaks_V3.pdf" ]]; then
        warn "diff file $f has a run tag; the email names the untagged files, edit it"
    fi
done

#
# Before this layout annotateWithHomer.sh wrote annote/ and
# R/diffAnalysisPairwise.R wrote its files in the analysis directory;
# such files are not delivered.
#
if [ -d annote ]; then
    warn "annote/ in the analysis directory is not delivered; move its files to $OUT/annote"
fi
shopt -s nullglob
OLD_DIFF=("${PROJ}"_*_DiffPeaks*_V3.*)
shopt -u nullglob
if [[ ${#OLD_DIFF[@]} -gt 0 ]]; then
    warn "${OLD_DIFF[*]} in the analysis directory are not delivered; move them to $DIFF"
fi

##
# Run folder: a trailing r_NNN as given, else the next one in ODIR

ODIR=$(realpath -m "$ODIR")
if [[ "$ODIR" =~ /(r_[0-9]+)$ ]]; then
    CURRDIR=${BASH_REMATCH[1]}
    ODIR=${ODIR%/*}
else
    CURRDIR=""
fi

#
# The deepest existing folder of ODIR must be readable by us; else
# the run folders in it could exist unseen, and r_001 would be picked
# or files added to a run folder without the note below.
#
seen=$ODIR
while [[ -n "$seen" && ! -e "$seen" ]]; do
    seen=${seen%/*}
done
seen=${seen:-/}
if [[ ! -r "$seen" || ! -x "$seen" ]]; then
    fatal "cannot read $seen, so cannot tell what is in $ODIR"
fi

if [[ -z "$CURRDIR" ]]; then
    last=0
    for d in "$ODIR"/r_[0-9]*; do
        n=${d##*/r_}
        if [[ -d "$d" && "$n" =~ ^[0-9]+$ ]] && ((10#$n > last)); then
            last=$((10#$n))
        fi
    done
    CURRDIR=$(printf 'r_%03d' $((last + 1)))
fi
DEST=$ODIR/$CURRDIR

echo "PROJECT=$PROJNO"
echo "GENOME=$GENOME"
echo "VERSION=$VERSION"
echo "ODIR=$ODIR"
echo "CURRDIR=$CURRDIR"
if [[ -e "$DEST" ]]; then
    echo "NOTE: $DEST exists; files are added to it"
fi
echo

##
# project.yaml first: a missing or unreadable README.txt or REQUEST.txt
# stops the delivery before anything is copied.

Rscript --no-init-file "$SDIR/bin/readme2yaml.R" "$PROJNO" "$ODIR" "$CURRDIR" "$GENOME"
echo

#
# The database takes pi and invest from README.txt; the delivery folder
# is usually named after the same two users.
#
YAML_PI=$(sed -n 's/^pi: //p' project.yaml)
YAML_INVEST=$(sed -n 's/^invest: //p' project.yaml)
if [[ "$ODIR" != "$DROOT/$YAML_PI/$YAML_INVEST/$PROJ" ]]; then
    warn "$ODIR is not $DROOT/<pi>/<invest>/$PROJ for pi=$YAML_PI invest=$YAML_INVEST of the request file"
fi

##
# Copy

asbot mkdir -p "$DEST/chipSeq/bw" "$DEST/chipSeq/docs"
asbot "${RSYNC[@]}" "${BW[@]}" "$DEST/chipSeq/bw/"
asbot "${RSYNC[@]}" "${FOLDERS[@]}" "$DEST/chipSeq/"
asbot "${RSYNC[@]}" "${DOCS[@]}" "$DEST/chipSeq/docs/"

##
# Email

EMAIL=deliveryEmail_${PROJNO}_${CURRDIR}_$(date +%y%m%d).txt

#
# A value as the replacement of sed s/// (git tags may contain / and &)
#
sedValue() {
    local v=${1//\\/\\\\}
    v=${v//\//\\/}
    printf '%s' "${v//&/\\&}"
}

fillTemplate() {
    sed -e "s/{PROJNO}/$(sedValue "$PROJNO")/g" \
        -e "s/{RUNDIR}/$(sedValue "$CURRDIR")/g" \
        -e "s/{VERSION}/$(sedValue "$VERSION")/g" \
        "$SDIR/assets/$1"
}

echo
echo "========================================================================="
echo
{
    fillTemplate tmpldeliveryEmail_01.txt
    if [[ "$HAVE_PEAKS" -eq 1 ]]; then
        fillTemplate tmpldeliveryEmail_02.txt
    fi
    if [[ "$HAVE_ANNOTE" -eq 1 ]]; then
        fillTemplate tmpldeliveryEmail_03.txt
    fi
    if [[ ${#DIFF_FILES[@]} -gt 0 ]]; then
        fillTemplate tmpldeliveryEmail_04.txt
    fi
    fillTemplate tmpldeliveryEmail_99.txt
} | tee "$EMAIL"
echo "========================================================================="
echo

##
# Permissions database

PERMS=("$BIC_DELIVERY/authorization_db/init_impact_project_permissions.py" -p project.yaml)

#
# Run with set +eu: the module scripts reference unset variables, and
# a failure here must say what was already done.
#
registerRun() {
    set +eu
    if [ "$(type -t module)" != "function" ]; then
        . /etc/profile.d/modules.sh || return 1
    fi
    module purge || return 1
    module load python/3.8.0 || return 1
    local python38
    python38=$(which python3.8) || return 1
    module load py-python-ldap/3.4.2 || return 1
    "$python38" "${PERMS[@]}"
}

if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRYRUN: module load python/3.8.0 py-python-ldap/3.4.2"
    echo "DRYRUN: python3.8 ${PERMS[*]@Q}"
else
    rc=0
    registerRun || rc=$?
    set -eu
    if [[ $rc -ne 0 ]]; then
        echo
        echo "    FATAL ERROR: the files are delivered to $DEST and the email"
        echo "    is in $EMAIL, but registering the run with the permissions"
        echo "    database failed (exit $rc). Fix the cause and run, here:"
        echo
        echo "        module load python/3.8.0 py-python-ldap/3.4.2"
        echo "        python3.8 ${PERMS[*]@Q}"
        echo
        exit 1
    fi
fi

echo
if [[ "$DRYRUN" -eq 1 ]]; then
    echo "DRY RUN: would deliver to $DEST; nothing was written there or to the database"
else
    echo "Delivered to $DEST"
fi
echo "Email: $EMAIL"
for w in "${WARNINGS[@]}"; do
    echo "    WARNING: $w"
done
echo
