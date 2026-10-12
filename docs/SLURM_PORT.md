# Porting the ChIP-seq pipeline from JUNO/LSF to IRIS/Slurm

Done 2026-10-11 on `feat/iris-port`, for v1.0.0. The design and its
decisions (D1 to D20) are in `00.SPEC_v1.0.0.md`. The model is the
ATAC-seq port (`../ATAC-seq/docs/SLURM_PORT.md`, v1.5.1 plus `devs`
28978c9): this pipeline had the same shape (a synchronous control job,
`xargs -n N bsub`, `bSync`), so the ATAC code transferred almost line for
line. JUNO is gone; there is no LSF path left (`attic/lsf.sh`,
`attic/bSync` are reference only).

## Strategy: keep the control job, swap the three primitives

| LSF | Slurm | Where |
| --- | --- | --- |
| `bsub -o LSF.CTRL/ ./pipe.sh` | `sbatch ./ChIP-seq/pipe.sh` | `~/bin/sbatch`, which exports `SBATCH_SCRIPT_DIR` |
| `bsub` | `~/bin/bsub` | SchedMD openlava shim with local patches |
| `bSync NAME` | `bSync.sh "^NAME$"` | `~/bin/bSync.sh`, wrapped in `bin/slurmTools.sh` |
| none | `bCheck NAME` (`sacct` State) | `bin/slurmTools.sh` |
| none | `bCheckAll REGEX` | `bin/slurmTools.sh` |

## Flag translation in pipe.sh

| LSF (v0.9.0) | Slurm (through the shim) | Note |
| --- | --- | --- |
| POST `-W 24:00 -R "rusage[mem=24]"` | `-n 2 -M 32G`, class from `runClass` | LSF `-W 24:00` is 24 hours; to `sbatch`/the shim `24:00` would be 24 minutes. 32G gives picard a 23g heap (`bin/picardJvm.sh`); 2 CPUs for its GC threads |
| BW `-n 4 -R "rusage[mem=10]"` (40G) | `-n 2 -M 24G` | `sort -S20g` plus predictd |
| CALLP `-W 359 -n 3 -R "rusage[mem=10]"` | `-n 3 -M 24G` | |
| MergePeaks `-n 3 -R "rusage[mem=10]"` | `-n 3 -M 24G`, `SHORT` | `sort -S 16g -T $TMPDIR` |
| Count `-R "rusage[mem=24]" -w "post_done(...)"`, `-T 10` | `-n 8 -M 24G`, threads from the job | `-w` removed: Slurm needs numeric ids, and an `afterok` on a failed parent pends forever. Stages 4 to 6 now run in sequence with `bSync` |
| none | CountCheck `-n 2 -M 8G` | Check B, one job per sample |
| DESEQ `-R "rusage[mem=24]"` | `-M 24G`, `SHORT` | |

`-M` is the job's total memory and a hard cgroup cap; LSF `rusage` was per
slot and advisory.

## Walltime classes

As ATAC: `SHORT` (`-W 1:59:00`, `cmobic_short,cpushort`, no qos) or `LONG`
(`-W 12:00:00`, `cmobic_cpu`, qos priority on that one `bsub` call). Stages
whose run time grows with input pick per job with `runClass MIN_PER_GB
FILE...`; `SHORT` if the estimate is at most 60 min
(`CHIP_SHORT_MAX_MIN`).

The rates are not yet measured for ChIP-seq. They are ATAC's, raised for
POST, and are to be replaced from the first full-size run (spec section 6,
item 11):

| Rate | Initial min/GB | Input | Reasoning |
| --- | --- | --- | --- |
| `RATE_POST` | 12 | input BAM | ATAC measured 3.5 to 5.3 for one picard sort; ChIP does two plus fixmate |
| `RATE_BW` | 30 | `clean.bed.gz` | ATAC's value |
| `RATE_CALLP` | 10 | target + control BED | ATAC's value |
| `RATE_COUNT` | 1.2 | all postProcess BAMs | ATAC's featureCounts 1.6.0 value; 2.1.1 is faster |
| `RATE_CCHECK` | 1.5 | postProcess BAM | ATAC's value |

## Run records and checking a run

As ATAC, with `CHIP` in place of `ATAC`: `00.RUNSTATUS.txt` (adds
`GENOME`, `PEAK_TYPE`, `PROPER_PAIR`, `SE`, `PAIRS`, `FRAGLEN`),
`SLURM.CTRL/jobs.tsv`, `#CHIP_EXIT=<rc>` trailers, `bin/checkRun.sh` (exit
0 ok, 1 failed, 2 running, 3 sacct failed). `MESSAGE` holds the reason for
a failure: the failed up-front check, or the first FATAL line of each
failed stage job (`failedJobsMessage`). `STATUS=COMPLETED` is written
last, after the `sacct` sweep and the deliverables check (spec section
3.3). One run per analysis directory; `00_PIPE_DONE` of LSF-era runs also
blocks a new run.

```bash
./ChIP-seq/bin/checkRun.sh
grep -L "#CHIP_EXIT=0" SLURM.0*/*.out
```

## Environment changes

| What | Was | Now |
| --- | --- | --- |
| samtools | `module load samtools` (unversioned; missing in the NoPP and SE scripts) | `module load samtools/1.20` in `bin/loadTools.sh`, sourced by every stage script |
| bedtools | `module load bedtools` | `~/bin/bedtools` 2.31.1; no module exists |
| picard | `bin/picardV2` loading the first module matching "17" (on IRIS `java/1.7.0_80`), `-Xmx23g`, `TMP_DIR=/scratch/socci` | `bin/picardV2`: `java/23.0.1`, heap from `SLURM_MEM_PER_NODE` minus 9G (`bin/picardJvm.sh`, from PEMapper), `TMP_DIR=$TMPDIR`; still picard 2.25.5 |
| POST scratch | relative `_scratch_ChIPSeq/<uuid>`, never removed | `${CHIP_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ChIPSeq}/<uuid>`, removed on exit |
| TMPDIR | node `/tmp`; MACS `/fscratch/socci` | per job under `${CHIP_LOCAL_TMP:-/localscratch/$USER}` (`bin/runStage.sh`) |
| venv | python 3.10.2 on `/juno` | `python3` on PATH (3.10.14), `numpy==1.23.0`, `MACS2==2.2.9.1` (`00.SETUP.sh`) |
| R | 4.2.2, user `~/.Rprofile` | 4.5.1, `Rscript --no-init-file` |
| HOMER | `/juno/bic/work/socci/opt/...`, hg19 only | `/usersoftware/core001/common/RHEL_8/homer`, hg19, hg38, mm10 |
| featureCounts | 2.0.6, `-p` only (counted reads) | 2.1.1, `-p --countReadPairs`, Checks A and B |

## Validation

Done 2026-10-11 on `feat/iris-port` (uncommitted, so the runs report
`v0.9.0-1-ge3d7535-dirty`), against `00.SPEC_v1.0.0.md` section 6. Test
directories, inputs and scripts are in `PortATACseq/PortChIPseq/`. The
downsampled set is the 8 Proj_16852 mm10 Cut&Run BAMs at 1/20
(`Ports/PEMapper/bin/downsampleBam.sh -s 20`, 0.23 to 0.35 GB each) in
`inputs/out/picard/DSB/`; GRCh38 BAMs are ATAC-seq's
`PortATACseq/out/picard/DSB/*/*.dn_100.bam`.

| Item | Test | Result |
| --- | --- | --- |
| 1 normal run | test01 (control job 18942608) | `COMPLETED`, 33 of 33 jobs `COMPLETED`, every log `#CHIP_EXIT=0`, `checkRun.sh` 0; Check A ratio 1.0000 and Check B diff 0 for all 8; FRAGLEN 147; 20 min wall, mostly queue |
| 2 matrix | `scripts/countTests.sh` on test01 | 2.0.6 with `-p --countReadPairs` gives the 2.1.1 matrix exactly; 2.0.6 with `-p` alone gives ratio 2.0000 and `checkCounts.sh` fails all 8. On the full-size Proj_16852 baseline summary `checkCounts.sh` also fails all 8 at 2.0000, finding the BAMs in `out/bam/` |
| 3 trap | test09, test11: copy of the checkout without `--countReadPairs` | `FAILED`, `STAGE=Count`, rc 1; test11 `MESSAGE=Count failed: ... FATAL ERROR: Check A failed: CHECK_A Cyc_1 ... ratio=2.0000 FAIL` |
| 4 no pairing file | test02 | `COMPLETED`, `MESSAGE=no --pairing-file; stages 3-6 skipped`, stage 1-2 deliverables present |
| 5 narrow, NoPP | test03 (`-n --proper-pair-off`) | `COMPLETED`, 33 jobs; `_peaks.narrowPeak`; Check A 1.0043 to 1.0047 (singletons kept) |
| 6 single-end | test04, fixture of 4 BAMs (read 1 only, `scripts/makeSEFixture.sh`) | `COMPLETED`, 17 jobs; featureCounts without `-p`; Check A in reads, 1.0000; Check B diff 0 |
| 7 failure injection | `scripts/failChecks.sh`; test05, test06, test10, test12; a checkout copy with a changed chr1 length | Before submission, rc 1, no jobs, reason in `MESSAGE`: mixed genomes, unknown genome, two SM tags, duplicate SM, paired BAMs with `-s`, single-end without `-s`, mixed read types, missing pairing file, missing BAM, header contig length differs. Usage errors (`-o`, bad option, no BAM) exit 1 without a status file (`00.ISSUES.md` #8). A second run in the same directory is refused and changes nothing. Pairing target `Cyc` (prefix of `Cyc_1`) and control `IgG_9` fail at `STAGE=CALLP` with the file line. A BAM on which predictd builds no model fails at `STAGE=BW` (test10). `scancel` of the control job: `FAILED`, rc 143, stage jobs `CANCELLED` (test12) |
| 8 repeat | `scripts/countTests.sh`, 20 runs at `-T 8` | 0 of 20 differ |
| 9 postPipe | test01 with Proj_16852 `manifest.txt` | QC xlsx and PDF in `out/qc/`, `___sigPeaks` and `___Volcano` in `out/qc/peaks/`; titles "Fragments in MACS Peaks", "Fraction of Fragments in MACS Peaks" |
| 10 Check B | `checkCounts.sh` on test01 with one result +1, one removed | rc 1, sample named; rc 1 on the count of results |
| 15 HOMER | test01 | 218,919 peaks, 99.9% with a gene name, no `chrchr`. Exits 1 with no status file, and on a peak file of the hg38 run (and the reverse): peaks beyond the contig lengths |
| 16 hg38 | test07, 4 ATAC BAMs, `na` controls | `COMPLETED`, `GENOME=hg38`; HOMER hg38 annotates 99.9% of 17,658 peaks; `diffAnalysisPairwise.R` runs (hg38 TxDb, 0 significant on this ATAC data) |
| final code | test13, as test01, after the last changes | `COMPLETED`, 33 of 33 jobs, 16 checks OK, count matrix identical to test01 |

Also run: `diffAnalysisPairwise.R` on test01 (mm10, 2 comparisons; the
significant peaks are annotated by ChIPseeker); it exits 1 for a manifest
`MapID` that is not a sample and outside a run directory.

Found during validation and fixed: `macs2 predictd` exits 0 when it cannot
build its model, so `makeBigWigFromBEDZ.sh` also requires an estimate in
the log (test08 had stopped only at `STAGE=FRAGLEN`).

Not done yet (spec section 6): item 11 to 14, the full-size Proj_16852
rerun (rates, comparison with the 2025 baseline); item 17, b37
(`00.ISSUES.md` #1).
