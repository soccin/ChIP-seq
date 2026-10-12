#
# JVM settings for bin/picardV2. Sourced, not executed. Copied from
# PEMapper v_4.1.0 bin/picardJvm.sh with the PEMAP_ variables renamed
# CHIP_.
#
# The old picardV2 ran a hardcoded -Xmx23g that had no connection to the
# memory the stage job asked for. Deriving the heap from the allocation
# keeps the two tied.
#

#
# Slurm exports SLURM_MEM_PER_NODE (in MB) whenever the job was submitted
# with --mem (bsub -M), which every POST job is. That number is a hard
# cgroup cap, so the heap has to sit below it by enough to cover what the
# JVM allocates outside the heap -- metaspace, thread stacks, direct
# buffers, GC bookkeeping -- and the page cache picard's own writes push
# into the same cgroup.
#
# The margin is fixed rather than proportional so that 32G gives back the
# historical -Xmx23g exactly: 32768 - 9216 = 23552. Unset, as when
# picardV2 is run by hand, falls back to that same 23g. Java's own cgroup
# awareness is no substitute: left to itself it takes 25% of the cap.
#
CHIP_XMX_MB=23552
CHIP_JVM_HEADROOM_MB=${CHIP_JVM_HEADROOM_MB:-9216}

if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    CHIP_XMX_MB=$(( SLURM_MEM_PER_NODE - CHIP_JVM_HEADROOM_MB ))
fi

if [ "$CHIP_XMX_MB" -lt 2048 ]; then
    echo "FATAL ERROR: derived picard heap [${CHIP_XMX_MB}m] is too small"
    echo "SLURM_MEM_PER_NODE=[${SLURM_MEM_PER_NODE:-}] headroom=[${CHIP_JVM_HEADROOM_MB}m]"
    echo "Raise the memory of the stage job, or lower CHIP_JVM_HEADROOM_MB."
    exit 1
fi

#
# A hard ceiling on GC threads, not a tuning knob: several picard JVMs on
# one node with their GCs sized to the node's physical core count will
# bury that node. POST asks for -n 2.
#
CHIP_JVM_GC_OPTS="-XX:ParallelGCThreads=2 -XX:ConcGCThreads=2"
