#!/bin/bash

#
# Control job for IRIS/Slurm. Run it from the analysis directory; every
# stage is fanned out with ~/bin/bsub, onto the short partitions unless
# its input is too large (see runClass), and the control job blocks
# between stages with bSync/bCheck. The control job itself can run for
# hours on full-size BAMs, so the directives below put it on the long
# partition with --qos=priority. Flags given to sbatch on the command
# line override them.
#
# CMD:
#    sbatch ./ChIP-seq/pipe.sh --pairing-file PAIRS BAM1 [BAM2 ...]
#
# Then, at any time, from the same directory:
#
#    ./ChIP-seq/bin/checkRun.sh
#
# which exits 0 if the run worked, 1 if it failed, 2 if it is still
# running and 3 if sacct failed and the state cannot be read now.
#
# One run per analysis directory. A run started in a directory that
# already holds one (see checkPreviousRun) stops at once and changes
# nothing; rerun in a new directory.
#
# sbatch here is the ~/bin wrapper. Slurm runs a copy of the batch script
# from its spool directory, so $0 cannot locate this checkout; the wrapper
# exports SBATCH_SCRIPT_DIR and SDIR is taken from that. Running pipe.sh
# directly (no sbatch) falls back to $0.
#
# Ported from the ATAC-seq control job (devs 28978c9); see
# docs/SLURM_PORT.md.
#
#SBATCH -p cmobic_cpu
#SBATCH --qos=priority
#SBATCH -t 3-00:00:00
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -J CTRL.ChIP
#SBATCH -o SLURM.CTRL/%j.out

set -eu

SDIR=${SBATCH_SCRIPT_DIR:-$( cd "$( dirname "$0" )" && pwd )}

if [ ! -e "$SDIR/postMapBamProcessing_ChIPSeq.sh" ]; then
    echo
    echo "    FATAL ERROR: cannot resolve SDIR=[$SDIR]"
    echo "    submit with ~/bin/sbatch, which exports SBATCH_SCRIPT_DIR"
    echo
    exit 1
fi

source $SDIR/bin/loadTools.sh
source $SDIR/bin/slurmTools.sh

if [ ! -e "$SDIR/venv" ]; then
    echo
    echo "   Need to install macs2"
    echo "   Source $SDIR/00.SETUP.sh in $SDIR (see README.md)"
    echo
    exit 1
fi

SCRIPT_VERSION=$(git --git-dir=$SDIR/.git --work-tree=$SDIR describe --always --long --tags --dirty)
PIPENAME="ChIP-Seq"

ORIG_CMD="$SDIR/pipe.sh $(printf '%q ' "$@")"

export PATH=$SDIR/bin:$PATH

##
# Process command args

function usage {
    cat <<-END_USAGE

	usage: sbatch ChIP-seq/pipe.sh [OPTIONS] --pairing-file <PAIRS> BAM1 [BAM2 ... BAMN]
	version=$SCRIPT_VERSION

	Options:

	  --pairing-file <PAIRS>   Tab/space delimited file of CONTROL/TARGET sample
	                           pairs. Required for peak calling (MACS); see below.

	  -n|--narrow-peaks        Call narrow peaks. Default is broad peaks.

	  --proper-pair-off        Keep reads that are not in a proper pair. Use for
	                           samples where translocations are important.

	  -s|--single-end-on       Input BAMs are single-end. Checked against the
	                           BAMs: required for single-end BAMs, refused for
	                           paired-end ones.

	Samples are named by the @RG SM tag of each BAM: one SM per BAM, unique
	across BAMs. Outputs are out/<SM>_postProcess.*. The genome is detected
	from the BAM headers; supported: $($SDIR/bin/genomeInfo.sh --tags).

	Pairing file:

	  Two columns, control in col-1 and target in col-2, one pair per line,
	  each a sample name (SM tag) exactly:

	      CTRL_1  TRGT_1
	      CTRL_2  TRGT_2

	  For target-only samples (no control) put 'na' in col-1:

	      na      TRGT_1
	      na      TRGT_2

	  WITHOUT --pairing-file the pipeline stops after the BAM post-processing
	  and bigWig stages; MACS peak calling and all downstream stages will NOT
	  be run.

	Example:

	  sbatch ./ChIP-seq/pipe.sh \\
	      --pairing-file pipeline/Proj_10706_B_sample_pairing.txt \\
	      pipeline/alignments/Proj_10706_B_s_*bam

	END_USAGE
    exit 1
}

PROPER_PAIR="Yes"
SE="No"
PAIRS=""
PEAK_OPT=""
PEAK_TYPE="broad"

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--narrow-peaks)
            PEAK_OPT="-n"
            PEAK_TYPE="narrow"
            shift
            ;;
        --proper-pair-off)
            PROPER_PAIR="No"
            shift
            ;;
        -s|--single-end-on)
            SE="Yes"
            shift
            ;;
        --pairing-file)
            if [[ $# -lt 2 ]]; then
                echo
                echo "    FATAL ERROR: --pairing-file needs a file"
                usage
            fi
            PAIRS="$2"
            shift 2
            ;;
        -o|--outdir)
            echo
            echo "    FATAL ERROR: -o|--outdir was removed in v1.0.0; output is always out/"
            usage
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo
            echo "    FATAL ERROR: invalid option [$1]"
            usage
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL_ARGS[@]}"

if [ "$#" -lt "1" ]; then
    usage
fi

TAG=q$PIPENAME

BAMS=$*

if [ "$PAIRS" == "" ]; then
    cat <<-END_WARNING

	#######################################################################

	  WARNING: no --pairing-file given

	  MACS peak calling will NOT be run. The pipeline will stop after
	  Stage2 (BAM post-processing and bigWig files); Stages 3-6 (peak
	  calling, peak merging, counts, scale factors) will be skipped.

	  To call peaks, rerun with --pairing-file <PAIRS>. For target-only
	  samples put 'na' in col-1 of the pairing file.

	#######################################################################

	END_WARNING
fi

echo SDIR=$SDIR
echo BAMS=$BAMS
echo VERSION=$SCRIPT_VERSION
echo PROPER_PAIR=$PROPER_PAIR SE=$SE PEAK_TYPE=$PEAK_TYPE PAIRS=$PAIRS

#
# Every sacct query in bCheck is bounded by this timestamp. Job names
# embed $$ and pids recycle, so without it a name from an earlier run
# could be counted against this one.
#
export CHIP_START=$(date -d '-2 min' +%Y-%m-%dT%H:%M:%S)
RUNRE="^${TAG}_.*_$$\$"

#
# Run records, read by bin/checkRun.sh together with sacct.
#
#   00.RUNSTATUS.txt     KEY=VALUE. STATUS is RUNNING from here on. It is
#                        set to COMPLETED only as the last action of this
#                        script, after the deliverables check, and to
#                        FAILED (with the STAGE it failed in) by the EXIT
#                        trap. If the control job is killed outright
#                        (SIGKILL on its memory cap, node failure) the
#                        file stays at RUNNING and checkRun.sh finds the
#                        control job dead in sacct. GENOME= is read by
#                        annotateWithHomer.sh and diffAnalysisPairwise.R.
#
#   SLURM.CTRL/jobs.tsv  one line per stage job: job id, stage, sample
#                        and log path, written at submission.
#
#   out/pipeline_info/   samples.tsv (SAMPLE, BAM, PBAM), read by every
#                        later step; cmd.sh.log, the provenance record.
#
RUNSTATUS=00.RUNSTATUS.txt
JOBS=SLURM.CTRL/jobs.tsv
INFODIR=out/pipeline_info
SAMPLES_TSV=$INFODIR/samples.tsv

CTRL_JOBID=${SLURM_JOB_ID:-none}
CTRL_LOG=none
if [ "$CTRL_JOBID" != "none" ]; then
    CTRL_LOG=$(scontrol show job $CTRL_JOBID 2>/dev/null | sed -n 's/^ *StdOut=//p') || true
    CTRL_LOG=${CTRL_LOG:-SLURM.CTRL/$CTRL_JOBID.out}
fi
echo CTRL_JOBID=$CTRL_JOBID

RUN_STATUS=RUNNING
RUN_STAGE=SETUP
RUN_RC=""
RUN_MESSAGE=""
RUN_STARTED=$(date '+%Y-%m-%d %H:%M:%S')
RUN_FINISHED=""
RUN_DONE=0
GENOME=""
SAMPLES=""
FRAGLEN=""
BSYNC_PID=""

writeRunStatus() {
    cat >$RUNSTATUS.tmp <<EOF
STATUS=$RUN_STATUS
STAGE=$RUN_STAGE
RC=$RUN_RC
MESSAGE=$RUN_MESSAGE
STARTED=$RUN_STARTED
FINISHED=$RUN_FINISHED
VERSION=$SCRIPT_VERSION
CTRL_JOBID=$CTRL_JOBID
CTRL_LOG=$CTRL_LOG
HOST=$(hostname)
PID=$$
CHIP_START=$CHIP_START
GENOME=$GENOME
PEAK_TYPE=$PEAK_TYPE
PROPER_PAIR=$PROPER_PAIR
SE=$SE
PAIRS=$PAIRS
FRAGLEN=$FRAGLEN
SAMPLES=$SAMPLES
BAMS=$BAMS
JOBS=$JOBS
EOF
    mv $RUNSTATUS.tmp $RUNSTATUS
}

setStage() {
    RUN_STAGE=$1
    writeRunStatus
    echo
    echo "==== STAGE $RUN_STAGE $(date '+%Y-%m-%d %H:%M:%S')"
}

#
# fatal MESSAGE
#
# Stop the run with MESSAGE in the status file (through the EXIT trap).
#
fatal() {
    RUN_MESSAGE="$*"
    echo
    echo "    FATAL ERROR: $RUN_MESSAGE"
    echo
    exit 1
}

#
# One run per analysis directory; pipe.sh is not built to be rerun. A
# second run would overwrite the first run's status file and job
# manifest, and the stages write into out/, so the earlier run's output
# would be mixed into this one. Stop before anything is written if the
# directory holds any of these. 00_PIPE_DONE is the stage sentinel of
# LSF-era runs. SLURM.CTRL itself is not a sign: Slurm creates it for
# this job's log. A control job that Slurm requeues after a node failure
# also stops here, on the records of its own first attempt. The check
# runs before the EXIT trap is installed, so the earlier run's records
# are left as they are; the reason is in this job's log.
#
checkPreviousRun() {
    local found=""
    local f

    for f in $RUNSTATUS $JOBS out 00_PIPE_DONE; do
        if [ -e "$f" ]; then
            found="$found $f"
        fi
    done

    if [ -n "$found" ]; then
        echo
        echo "    FATAL: this directory already holds a run:$found"
        echo "    pipe.sh cannot be rerun in the same directory; use a new one."
        echo
        exit 1
    fi
}

checkPreviousRun

mkdir -p SLURM.CTRL
printf '#JOBID\tSTAGE\tSAMPLE\tLOG\n' >$JOBS
writeRunStatus

#
# Any exit before RUN_DONE=1 is a failure: an error under set -e, a failed
# bCheck, a failed check (fatal), or scancel / time limit (SIGTERM).
# Record it in the status file and cancel the rest of this run's stage
# jobs so they do not keep running after the control job is gone. SIGTERM
# reaches the trap promptly only while the control job is in bSync (see
# slurmTools.sh), which is where it spends nearly all of its time.
#
onSignal() {
    RUN_MESSAGE="control job got SIG$1 (scancel, time limit or node shutdown)"
    exit $2
}
trap 'onSignal TERM 143' TERM
trap 'onSignal INT 130' INT

onExit() {
    local rc=$?
    if [ -n "$BSYNC_PID" ]; then
        kill $BSYNC_PID 2>/dev/null || true
    fi
    if [ "$RUN_DONE" != "1" ]; then
        #
        # A failed bCheck leaves the reason in the stage job logs; copy it
        # into MESSAGE.
        #
        local reason=""
        reason=$(failedJobsMessage "$RUN_STAGE" "$JOBS" 2>/dev/null) || true
        if [ -n "$reason" ]; then
            RUN_MESSAGE="$reason${RUN_MESSAGE:+; $RUN_MESSAGE}"
        fi
        RUN_STATUS=FAILED
        RUN_RC=$rc
        RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
        writeRunStatus
        echo
        echo "    pipe.sh FAILED in stage $RUN_STAGE (rc=$rc); cancelling queued jobs of this run"
        echo
        bKill "$RUNRE"
    fi
}
trap onExit EXIT

##############################################################################
# Up-front checks. Nothing is submitted until all of them pass.
##############################################################################

setStage CHECK

if [ -n "$PAIRS" ] && [ ! -r "$PAIRS" ]; then
    fatal "cannot read pairing file [$PAIRS]"
fi

for BAM in "$@"; do
    if [ ! -r "$BAM" ]; then
        fatal "cannot read BAM [$BAM]"
    fi
done

#
# Genome. Every BAM is identified by bin/getGenomeBuildBAM.sh and all must
# agree; the build must be a row of lib/genomes/GENOMES.tsv; and every
# contig of lib/genomes/<GENOME>.genome must be in every BAM header with
# the same length. The last check catches a header md5 mapped to the
# wrong tag and a genome file that does not match the data.
#
SUPPORTED_GENOMES=$($SDIR/bin/genomeInfo.sh --tags)

#
# GENOME is set (and so recorded) only once every check has passed.
#
FIRST_BUILD=""
for BAM in "$@"; do
    BUILD=$($SDIR/bin/getGenomeBuildBAM.sh $BAM) \
        || fatal "cannot identify the genome of [$BAM]"
    echo "    genome [$BUILD] $BAM"
    if [ -z "$FIRST_BUILD" ]; then
        FIRST_BUILD=$BUILD
    elif [ "$BUILD" != "$FIRST_BUILD" ]; then
        fatal "BAMs on different genomes: [$FIRST_BUILD] ($1) and [$BUILD] ($BAM)"
    fi
done

case " $SUPPORTED_GENOMES " in
    *" $FIRST_BUILD "*)
        ;;
    *)
        fatal "unsupported genome [$FIRST_BUILD]; supported: $SUPPORTED_GENOMES"
        ;;
esac

GENOME_FILE=$($SDIR/bin/genomeInfo.sh $FIRST_BUILD GENOME_FILE) \
    || fatal "no genome file for [$FIRST_BUILD]"

for BAM in "$@"; do
    MISSING_CONTIGS=$(samtools view -H $BAM \
        | awk -F'\t' '
            NR == FNR {len[$1] = $2; next}
            $1 == "@SQ" {
                sn = ""; ln = ""
                for (i = 2; i <= NF; i++) {
                    if ($i ~ /^SN:/) sn = substr($i, 4)
                    if ($i ~ /^LN:/) ln = substr($i, 4)
                }
                hdr[sn] = ln
            }
            END {
                for (c in len) if (hdr[c] != len[c]) printf "%s ", c
            }' $GENOME_FILE -) \
        || fatal "cannot read the header of [$BAM]"
    if [ -n "$MISSING_CONTIGS" ]; then
        fatal "header of [$BAM] lacks contigs of $GENOME_FILE, or their lengths differ: $MISSING_CONTIGS"
    fi
done

GENOME=$FIRST_BUILD
echo GENOME=$GENOME
writeRunStatus

#
# Samples. The sample name is the @RG SM tag: exactly one distinct value
# per BAM, unique across BAMs. It names every output (out/<SM>_postProcess.*)
# and is what the pairing file and the manifests must name, exactly. It
# must be a plain file name, and 'na' means "no control" in the pairing
# file, so it cannot be a sample name.
#
getSMTags() {
    samtools view -H "$1" \
        | awk -F'\t' '$1 == "@RG" {
            for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4)
          }' \
        | sort -u
}

declare -A SAMPLE_OF
for BAM in "$@"; do
    SM=$(getSMTags $BAM)
    NSM=$(echo "$SM" | awk 'NF > 0' | wc -l)
    if [ "$NSM" == "0" ]; then
        fatal "no @RG SM tag in [$BAM]"
    fi
    if [ "$NSM" != "1" ]; then
        fatal "more than one SM tag in [$BAM]: $(echo $SM)"
    fi
    if [[ ! "$SM" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        fatal "SM tag [$SM] of [$BAM] is not a plain name (letters, digits, . _ -)"
    fi
    if [ "${SM,,}" == "na" ]; then
        fatal "SM tag [$SM] of [$BAM]: 'na' means no control in the pairing file"
    fi
    SAMPLE_OF[$BAM]=$SM
    SAMPLES="$SAMPLES $SM"
done
SAMPLES=${SAMPLES# }

DUPS=$(echo $SAMPLES | tr ' ' '\n' | sort | uniq -d)
if [ -n "$DUPS" ]; then
    fatal "SM tag used by more than one BAM: $(echo $DUPS)"
fi

echo SAMPLES=$SAMPLES
writeRunStatus

#
# Read type, from the first 10,000 records of each BAM: paired (flag 0x1
# on every record) or single-end (on none). All BAMs must agree, and -s
# must be given exactly when they are single-end; otherwise paired BAMs
# would go through the single-end recipe and be counted as reads.
#
readType() {
    samtools view "$1" 2>/dev/null \
        | head -n 10000 \
        | awk '{if (and($2, 1)) p++; else s++}
               END {
                   if (p + s == 0) print "empty"
                   else if (s == 0) print "paired"
                   else if (p == 0) print "single"
                   else print "mixed"
               }'
}

READ_TYPE=""
for BAM in "$@"; do
    RT=$(readType $BAM)
    echo "    reads [$RT] $BAM"
    case $RT in
        paired|single)
            ;;
        *)
            fatal "cannot tell paired from single-end reads in [$BAM]: $RT"
            ;;
    esac
    if [ -z "$READ_TYPE" ]; then
        READ_TYPE=$RT
    elif [ "$RT" != "$READ_TYPE" ]; then
        fatal "BAMs mix paired-end and single-end reads: [$BAM] is $RT, [$1] is $READ_TYPE"
    fi
done

if [ "$READ_TYPE" == "paired" ] && [ "$SE" == "Yes" ]; then
    fatal "-s|--single-end-on given but the BAMs are paired-end"
fi
if [ "$READ_TYPE" == "single" ] && [ "$SE" == "No" ]; then
    fatal "the BAMs are single-end; run with -s|--single-end-on"
fi

if [ "$SE" == "Yes" ]; then
    POST_SCRIPT=postMapBamProcessing_ChIPSeq_SE.sh
    COUNT_OPT="--single-end"
elif [ "$PROPER_PAIR" == "Yes" ]; then
    POST_SCRIPT=postMapBamProcessing_ChIPSeq.sh
    COUNT_OPT=""
else
    POST_SCRIPT=postMapBamProcessing_ChIPSeq_NoPP.sh
    COUNT_OPT=""
fi
echo POST_SCRIPT=$POST_SCRIPT

#
# The sample table every later step reads. BAM is absolute; PBAM is the
# postProcess BAM as featureCounts names its column.
#
mkdir -p out $INFODIR
printf 'SAMPLE\tBAM\tPBAM\n' >$SAMPLES_TSV
for BAM in "$@"; do
    ABAM=$(cd "$(dirname "$BAM")" && pwd)/$(basename "$BAM")
    printf '%s\t%s\t%s\n' "${SAMPLE_OF[$BAM]}" "$ABAM" "out/${SAMPLE_OF[$BAM]}_postProcess.bam" >>$SAMPLES_TSV
done

sampleBAM() {
    awk -F'\t' -v s="$1" 'NR > 1 && $1 == s {print $2}' $SAMPLES_TSV
}

#
# Provenance, written at the end of every run that gets this far,
# whether or not it has a pairing file.
#
writeCmdLog() {
    local CMD_LOG=$INFODIR/cmd.sh.log
    local GURL
    GURL=$(git --git-dir=$SDIR/.git --work-tree=$SDIR config --get remote.origin.url) || GURL=""

    cat <<-END_VERSION >$CMD_LOG
	DATE: $(date)
	SDIR: $SDIR
	GURL: $GURL
	GTAG: $SCRIPT_VERSION
	PWD: $PWD
	CTRL_JOBID: $CTRL_JOBID

	PROPER_PAIR: $PROPER_PAIR
	PEAK_TYPE: $PEAK_TYPE
	SE: $SE
	PAIRS: $PAIRS
	GENOME: $GENOME
	FRAGLEN: $FRAGLEN
	SAMPLES: $SAMPLES
	BAMS: $BAMS

	TOOLS:
	  $($SDIR/bin/featureCounts -v 2>&1 | grep -m1 featureCounts || true)
	  $(samtools --version | head -1 || true)
	  $(bedtools --version || true)
	  $($SDIR/venv/bin/macs2 --version 2>&1 || true)
	  $(R --version | head -1 || true)
	  picard $(basename $(dirname $SDIR/bin/jar/*/picard.jar))

	Script: $ORIG_CMD
	END_VERSION
}

##############################################################################
# Stages
##############################################################################

#
# Scrub the Slurm variables of the control job itself. bsub submits with
# the default --export=ALL, so they would otherwise ride along into every
# stage job. SLURM_CONF stays; the client tools need it. The SBATCH_*
# input variables go too (SDIR has already been taken from
# SBATCH_SCRIPT_DIR): one left in the submitting shell would apply to
# every stage job, and SBATCH_QOS=priority there gets every SHORT job
# rejected. Stage partition, walltime and qos come only from chipSub.
#
for V in ${!SLURM_@}; do
    if [ "$V" != "SLURM_CONF" ]; then
        unset $V
    fi
done
for V in ${!SBATCH_@}; do
    unset $V
done

#
# Walltime classes. Every stage job is SHORT or LONG.
#
#   SHORT  -W 1:59:00, no qos. ~/bin/bsub sends anything under two hours
#          to cmobic_short,cpushort. cpushort allows only qos=normal and
#          EnforcePartLimits=ALL is set, so a SHORT job must not carry
#          a qos or it is rejected ("Invalid qos specification").
#   LONG   -W 12:00:00 and qos priority. bsub sends it to cmobic_cpu;
#          SBATCH_QOS=priority is set on that one bsub call only.
#
# SHORT whenever possible: jobs on cmobic_cpu can wait hours to start.
# Stages whose run time grows with the size of their input pick their
# class per job with runClass; the others are always SHORT. -M is the
# total memory for the job and a hard cgroup limit.
#
SHORT_W=1:59:00
LONG_W=12:00:00

#
# Minutes per GB of input. Not yet measured for ChIP-seq: these are the
# ATAC-seq rates, raised for POST (two picard sorts plus fixmate), and
# are to be replaced by sacct Elapsed from the first full-size run. See
# docs/SLURM_PORT.md.
#
RATE_POST=12       # input BAM
RATE_BW=30         # clean.bed.gz
RATE_CALLP=10      # target + control clean.bed.gz
RATE_COUNT=1.2     # all postProcess BAMs, summed
RATE_CCHECK=1.5    # postProcess BAM

#
# A job is SHORT if its estimated run time is at most SHORT_MAX_MIN, half
# the SHORT walltime, to leave room for slower nodes and a loaded shared
# filesystem. CHIP_SHORT_MAX_MIN overrides it (0 makes every estimated
# stage LONG).
#
SHORT_MAX_MIN=${CHIP_SHORT_MAX_MIN:-60}

#
# runClass MIN_PER_GB FILE [FILE ...]
#
# Print SHORT or LONG for a job whose input is FILE..., and log the
# estimate on stderr.
#
runClass() {
    local rate=$1
    local bytes=0
    local f size
    shift

    for f in "$@"; do
        size=$(stat -L -c %s "$f") || return 1
        bytes=$((bytes + size))
    done

    awk -v b=$bytes -v r=$rate -v m=$SHORT_MAX_MIN 'BEGIN {
        est = b / 1e9 * r
        class = (est <= m) ? "SHORT" : "LONG"
        printf "    runClass %.1f GB x %s min/GB = %.0f min -> %s\n", \
            b / 1e9, r, est, class >"/dev/stderr"
        print class
    }'
}

#
# chipSub STAGE SAMPLE LOGDIR CLASS "BSUB_OPTS" CMD [ARG ...]
#
# Submit one stage job as ${TAG}_STAGE_$$ in walltime class CLASS (SHORT
# or LONG) through bin/runStage.sh, which appends #CHIP_EXIT=<rc> to the
# job log, and record it in $JOBS.
#
chipSub() {
    local stage=$1
    local sample=$2
    local logdir=$3
    local class=$4
    local opts=$5
    local out jobid runtime
    local qos=()
    shift 5

    case $class in
        SHORT)
            runtime=$SHORT_W
            ;;
        LONG)
            runtime=$LONG_W
            qos=(SBATCH_QOS=priority)
            ;;
        *)
            fatal "bad class [$class] for stage $stage sample $sample"
            ;;
    esac

    if ! out=$(env "${qos[@]}" bsub -W $runtime -o $logdir/ \
                   -J ${TAG}_${stage}_$$ $opts \
                   $SDIR/bin/runStage.sh "$@"); then
        echo "$out"
        fatal "bsub failed for stage $stage sample $sample"
    fi
    echo "$out"

    jobid=$(echo "$out" | sed -n 's/^Job <\([0-9][0-9]*\)> is submitted.*/\1/p')
    if [ -z "$jobid" ]; then
        fatal "no job id in bsub output for stage $stage sample $sample"
    fi

    printf '%s\t%s\t%s\t%s\n' "$jobid" "$stage" "$sample" "$logdir/$jobid.out" >>$JOBS
}

#
# waitStage STAGE
#
# Block until every ${TAG}_STAGE_$$ job is done, then abort the run if any
# of them did not reach COMPLETED. The number of jobs submitted, from
# $JOBS, lets bCheck wait for every accounting record.
#
waitStage() {
    local njobs
    setStage $1
    njobs=$(awk -F'\t' -v s="$1" '$2 == s' $JOBS | wc -l)
    bSync ${TAG}_$1_$$
    bCheck ${TAG}_$1_$$ $njobs
}

#
# Deliverables. A peak file may legitimately be empty, so it only has to
# exist; everything else has to be non-empty.
#
MISSING=""

checkFile() {
    if [ ! -s "$1" ]; then
        MISSING="$MISSING $1"
    fi
}

checkExists() {
    if [ ! -e "$1" ]; then
        MISSING="$MISSING $1"
    fi
}

checkSampleDeliverables() {
    local S
    for S in $SAMPLES; do
        checkFile out/${S}_postProcess.bam
        if [ ! -s out/${S}_postProcess.bai ] && [ ! -s out/${S}_postProcess.bam.bai ]; then
            MISSING="$MISSING out/${S}_postProcess.bai"
        fi
        checkFile out/${S}_postProcess.clean.bed.gz
        checkFile out/profiles/${S}_postProcess.clean.10mNorm.bw
        checkFile out/profiles/${S}_postProcess.clean.log
    done
    checkFile $SAMPLES_TSV
    checkFile $INFODIR/cmd.sh.log
}

failIfMissing() {
    if [ -n "$MISSING" ]; then
        echo
        echo "    ERROR missing deliverables:"
        for F in $MISSING; do
            echo "        $F"
        done
        echo
        fatal "missing deliverables:$MISSING"
    fi
}

#
# finishRun MESSAGE
#
# The sweep of every stage job of the run, then COMPLETED.
#
finishRun() {
    local failed
    setStage STAGING
    failed=$(bCheckAll "$RUNRE")
    if [ "$failed" != "" ]; then
        echo -e "\n\n\nFailed Slurm jobs\n\n"
        echo "$failed"
        echo -e "\n\n"
        fatal "stage jobs not COMPLETED:$(echo "$failed" | awk '{printf " %s %s", $1, $2}')"
    fi

    RUN_STATUS=COMPLETED
    RUN_STAGE=DONE
    RUN_RC=0
    RUN_MESSAGE="$1"
    RUN_FINISHED=$(date '+%Y-%m-%d %H:%M:%S')
    writeRunStatus
    RUN_DONE=1

    echo
    echo "==== ChIP-Seq run COMPLETED $RUN_FINISHED"
    echo
}

setStage SUBMIT

#
# Stage 1: POST, one job per BAM.
#
for S in $SAMPLES; do
    BAM=$(sampleBAM $S)
    CLASS=$(runClass $RATE_POST $BAM) || fatal "cannot size [$BAM]"
    chipSub POST $S SLURM.01.POST $CLASS "-n 2 -M 32G" \
        $SDIR/$POST_SCRIPT $GENOME $S $BAM
done

waitStage POST

#
# Stage 2: BW, one job per clean BED: bigWig and macs2 predictd.
#
for S in $SAMPLES; do
    BEDZ=out/${S}_postProcess.clean.bed.gz
    CLASS=$(runClass $RATE_BW $BEDZ) || fatal "no POST output [$BEDZ]"
    chipSub BW $S SLURM.02.BW $CLASS "-n 2 -M 24G" \
        $SDIR/makeBigWigFromBEDZ.sh $GENOME $BEDZ
done

waitStage BW

#
# Stage 2b: one fragment length for every MACS call, the median of the
# per-sample predictd estimates.
#
setStage FRAGLEN

PREDICTD_LOGS=""
for S in $SAMPLES; do
    PREDICTD_LOGS="$PREDICTD_LOGS out/profiles/${S}_postProcess.clean.log"
done

FRAGLEN=$(Rscript --no-init-file $SDIR/getMedianFragmentLengthFromPredictDFile.R $PREDICTD_LOGS) \
    || fatal "no fragment length estimate (getMedianFragmentLengthFromPredictDFile.R failed)"
FRAGLEN=$(echo $FRAGLEN)

if [[ ! "$FRAGLEN" =~ ^[1-9][0-9]*$ ]]; then
    fatal "no fragment length estimate [$FRAGLEN]"
fi

echo FRAGLEN=$FRAGLEN
writeRunStatus

if [ "$PAIRS" == "" ]; then
    echo
    echo "No --pairing-file given; skipping MACS peak calling"
    echo "Stopping here. Stages 3-6 will not be run"
    echo "To run unpaired samples put 'na' in col 1"
    echo "of pairing file"
    echo
    writeCmdLog
    setStage DELIVERABLES
    checkSampleDeliverables
    failIfMissing
    finishRun "no --pairing-file; stages 3-6 skipped"
    exit 0
fi

#
# Stage 3: CALLP, one job per line of out/macsPairs.txt. A name in the
# pairing file that is not a sample is fatal here, with the line.
#
setStage CALLP

GMA_ERR=$(mktemp -p "$TMPDIR" generateMACSArgs.XXXXXX)
if ! Rscript --no-init-file $SDIR/generateMACSArgs.R $PAIRS $SAMPLES_TSV \
        >out/macsPairs.txt 2>$GMA_ERR; then
    cat $GMA_ERR
    fatal "generateMACSArgs.R: $(grep -m1 'FATAL' $GMA_ERR | tr -s '[:space:]' ' ' | cut -c1-300)"
fi
cat $GMA_ERR
rm -f $GMA_ERR

cat out/macsPairs.txt

mapfile -t MACS_PAIRS <out/macsPairs.txt
if [ "${#MACS_PAIRS[@]}" == "0" ]; then
    fatal "no pairs in out/macsPairs.txt"
fi

for LINE in "${MACS_PAIRS[@]}"; do
    read -r TBED CBED <<<"$LINE"
    TNAME=$(basename $TBED _postProcess.clean.bed.gz)
    if [ "$CBED" == "_NoCTRL" ]; then
        CLASS=$(runClass $RATE_CALLP $TBED) || fatal "cannot size [$TBED]"
    else
        CLASS=$(runClass $RATE_CALLP $TBED $CBED) || fatal "cannot size [$TBED $CBED]"
    fi
    chipSub CALLP $TNAME SLURM.03.CALLP $CLASS "-n 3 -M 24G" \
        $SDIR/callPeaks_ChIPseq.sh $PEAK_OPT $GENOME $FRAGLEN $TBED $CBED
done

waitStage CALLP

#
# MergePeaks -> Count -> DESEQ run in sequence through bSync rather than
# with -w post_done(): the bsub shim maps that to afterok but Slurm needs
# numeric job ids, and an afterok whose parent failed pends forever in
# DependencyNeverSatisfied.
#
SAF=out/macs/macsPeaksMerged.saf
COUNTS=out/macs/peaks_raw_fcCounts.txt

chipSub MergePeaks all SLURM.04.MERGE SHORT "-n 3 -M 24G" \
    $SDIR/mergePeaksToSAF.sh out/macs \>$SAF

waitStage MergePeaks

if [ ! -s "$SAF" ]; then
    fatal "no merged peaks in $SAF"
fi

#
# The counts are never used unchecked. featureCounts 2.0.6, run here with
# -p but without --countReadPairs from v0.7.0 to v0.9.0, counted every
# read of every BAM and exited 0. Two checks must pass before DESEQ:
#
#   A  in the Count job (bin/countPeaks.sh): each sample's fragment total
#      is about half its mapped reads (equal to them for -s).
#   B  one CountCheck job per sample (bin/assignedCheck.sh), run beside
#      Count: the fragments over the merged peaks, counted with samtools
#      alone, must equal featureCounts' Assigned. Compared here by
#      bin/checkCounts.sh once both stages are done.
#
# The POST scripts write the .bai (picard CREATE_INDEX), so there is no
# Index stage.
#
PBAMS=""
ATLAS_FRAGS=""
for S in $SAMPLES; do
    PBAMS="$PBAMS out/${S}_postProcess.bam"
    ATLAS_FRAGS="$ATLAS_FRAGS out/${S}_postProcess.atlasFragments.txt"
done

CLASS=$(runClass $RATE_COUNT $PBAMS) || fatal "cannot size the postProcess BAMs"
chipSub Count all SLURM.05.COUNT $CLASS "-n 8 -M 24G" \
    $SDIR/bin/countPeaks.sh $COUNT_OPT 10 $SAF $COUNTS $PBAMS

for S in $SAMPLES; do
    PBAM=out/${S}_postProcess.bam
    CLASS=$(runClass $RATE_CCHECK $PBAM) || fatal "cannot size [$PBAM]"
    chipSub CountCheck $S SLURM.05b.CHECK $CLASS "-n 2 -M 8G" \
        $SDIR/bin/assignedCheck.sh $PBAM $SAF out/${S}_postProcess.atlasFragments.txt
done

waitStage Count
waitStage CountCheck

CHECK_REPORT=out/macs/peaks_raw_fcCounts.checks.txt
if ! CHECK_OUT=$($SDIR/bin/checkCounts.sh $COUNT_OPT $COUNTS.summary $ATLAS_FRAGS); then
    echo "$CHECK_OUT"
    echo "$CHECK_OUT" >$CHECK_REPORT
    fatal "count check failed: $(echo "$CHECK_OUT" | grep -m1 FAIL || true)"
fi
echo "$CHECK_OUT"
echo "$CHECK_OUT" >$CHECK_REPORT

#
# Stage 6: DESeq2 scale factors.
#
chipSub DESEQ all SLURM.06.DESEQ SHORT "-M 24G" \
    Rscript --no-init-file $SDIR/getDESeqScaleFactors.R $COUNTS

waitStage DESEQ

writeCmdLog

setStage DELIVERABLES

checkSampleDeliverables

for LINE in "${MACS_PAIRS[@]}"; do
    read -r TBED CBED <<<"$LINE"
    TNAME=$(basename $TBED _postProcess.clean.bed.gz)
    if [ "$CBED" == "_NoCTRL" ]; then
        CNAME=___NoCTRL
    else
        CNAME=$(basename $CBED _postProcess.clean.bed.gz)
    fi
    PREFIX=${TNAME}_${CNAME}
    checkFile out/macs/$PREFIX/${PREFIX}_peaks.xls
    if [ "$PEAK_TYPE" == "narrow" ]; then
        checkExists out/macs/$PREFIX/${PREFIX}_peaks.narrowPeak
    else
        checkExists out/macs/$PREFIX/${PREFIX}_peaks.broadPeak
    fi
done

checkFile $SAF
checkFile $COUNTS
checkFile $COUNTS.summary
checkFile $CHECK_REPORT
checkFile out/macs/scaleFactorsDESeq2.csv
checkFile out/macsPairs.txt

failIfMissing

finishRun ""

cat <<EOF

pipe.sh done

  Check the run with

      $SDIR/bin/checkRun.sh

  then run

      $SDIR/postPipe.sh

EOF
