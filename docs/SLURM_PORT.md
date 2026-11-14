# Porting the ChIP-seq pipeline from JUNO/LSF to IRIS/Slurm

Done 2026-10-11 on `feat/iris-port`, for v1.0.0. The design and its
decisions (D1 to D20) are in `attic/00.SPEC_v1.0.0.md`. The model is the
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
| `bSync NAME` | `bSync "^NAME$"` | `~/bin/bSync`, wrapped in `bin/slurmTools.sh` |
| none | `bCheck NAME` (`sacct` State) | `bin/slurmTools.sh` |
| none | `bCheckAll REGEX` | `bin/slurmTools.sh` |

## Flag translation in pipe.sh

| LSF (v0.9.0) | Slurm (through the shim) | Note |
| --- | --- | --- |
| POST `-W 24:00 -R "rusage[mem=24]"` | `-n 2 -M 40G`, class from `runClass` | LSF `-W 24:00` is 24 hours; to `sbatch`/the shim `24:00` would be 24 minutes. Picard keeps a 23g heap (`bin/picardJvm.sh`, 40G minus 17G); the other 17G is headroom, since at 32G all four full-size POST jobs reached MaxRSS 32.0 GiB. 2 CPUs for its GC threads |
| BW `-n 4 -R "rusage[mem=10]"` (40G) | `-n 2 -M 24G` | `sort -S20g` plus predictd |
| CALLP `-W 359 -n 3 -R "rusage[mem=10]"` | `-n 3 -M 24G` | |
| MergePeaks `-n 3 -R "rusage[mem=10]"` | `-n 3 -M 24G`, `SHORT` | `sort -S 16g -T $TMPDIR` |
| Count `-R "rusage[mem=24]" -w "post_done(...)"`, `-T 10` | `-n 8 -M 24G`, threads from the job | `-w` removed: Slurm needs numeric ids, and an `afterok` on a failed parent pends forever. Stages 4 to 6 now run in sequence with `bSync` |
| none | CountCheck `-n 2 -M 16G` | Check B, one job per sample. Was 8G; the full-size jobs reached MaxRSS 7.3 to 8.0 GiB |
| DESEQ `-R "rusage[mem=24]"` | `-M 24G`, `SHORT` | |

`-M` is the job's total memory and a hard cgroup cap; LSF `rusage` was per
slot and advisory.

## Walltime classes

As ATAC: `SHORT` (`-W 1:59:00`, `cmobic_short,cpushort`, no qos) or `LONG`
(`-W 12:00:00`, `cmobic_cpu`, qos priority on that one `bsub` call). Stages
whose run time grows with input pick per job with `runClass MIN_PER_GB
FILE...`; `SHORT` if the estimate is at most 60 min
(`CHIP_SHORT_MAX_MIN`).

The rates were first ATAC's (POST raised to 12) and were checked against
`sacct` Elapsed of the full-size Proj_15235_B run (2026-11-13, 4 mm10
BAMs of 4.3 to 5.2 GB; acceptance item 11). POST and CCHECK were too low
and were raised; the others were above what the run measured and were
kept.

| Rate | min/GB | Input | Measured (Proj_15235_B) | Note |
| --- | --- | --- | --- | --- |
| `RATE_POST` | 25 (was 12) | input BAM | 17.6 to 20.7 (75 to 100 min) | At 12, two jobs ran SHORT; one ended 19 min before the 1:59 limit |
| `RATE_BW` | 30 | `clean.bed.gz` | 18.9 to 19.9 (14 to 20 min) | ATAC's value |
| `RATE_CALLP` | 10 | target + control BED | about 5 (4 to 5 min) | ATAC's value |
| `RATE_COUNT` | 1.2 | all postProcess BAMs | 0.07 (1 min for 13.6 GB, 8 threads, MaxRSS 14 GiB) | ATAC's value; kept at 8 threads (D14) |
| `RATE_CCHECK` | 2.5 (was 1.5) | postProcess BAM | 1.75 to 1.97 (5.5 to 7 min) | |

MergePeaks (1 s) and DESEQ (10 s) are always `SHORT`.

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
| picard | `bin/picardV2` loading the first module matching "17" (on IRIS `java/1.7.0_80`), `-Xmx23g`, `TMP_DIR=/scratch/socci` | `bin/picardV2`: `java/23.0.1`, heap from `SLURM_MEM_PER_NODE` minus 17G (`bin/picardJvm.sh`, from PEMapper, which uses 9G), `TMP_DIR=$TMPDIR`; still picard 2.25.5 |
| POST scratch | relative `_scratch_ChIPSeq/<uuid>`, never removed | `${CHIP_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ChIPSeq}/<uuid>`, removed on exit |
| TMPDIR | node `/tmp`; MACS `/fscratch/socci` | per job under `${CHIP_LOCAL_TMP:-/localscratch/$USER}` (`bin/runStage.sh`) |
| venv | python 3.10.2 on `/juno` | `python3` on PATH (3.10.14), `numpy==1.23.0`, `MACS2==2.2.9.1` (`00.SETUP.sh`) |
| R | 4.2.2, user `~/.Rprofile` | 4.5.1, `Rscript --no-init-file` |
| HOMER | `/juno/bic/work/socci/opt/...`, hg19 only | `/usersoftware/core001/common/RHEL_8/homer`, hg19, hg38, mm10 |
| featureCounts | 2.0.6, `-p` only (counted reads) | 2.1.1, `-p --countReadPairs`, Checks A and B |

## Validation

Done 2026-10-11 on `feat/iris-port` (uncommitted, so the runs report
`v0.9.0-1-ge3d7535-dirty`), against `attic/00.SPEC_v1.0.0.md` section 6. Test
directories, inputs and scripts were in `PortATACseq/PortChIPseq/`, since
removed. The
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

### Full-size run: Proj_15235_B (items 11 to 14)

Done 2026-11-13 on `chore/resume` (`v0.9.0-27-geffcb97`), in
`Proj_15235_B/ChIP_IRIS/`: the 4 `Map/out___M/*/*___MD.bam` (19.3 GB),
the baseline's `pairing.txt` (`na` controls), default options (broad,
proper pairs). Baseline: `Proj_15235_B/ChIP/out`, v0.9.0 on JUNO,
2026-07-19. It replaces the Proj_16852 rerun planned at first; a fuller
test on Proj_16852 is `00.ISSUES.md` #15. The comparison scripts are in
`ChIP_IRIS/validate/` (`item12_upstream.sh`, `item13_recount.sh`,
`compareCounts.R`, `item14_postPipe.sh`), with their reports beside
them.

| Item | Result |
| --- | --- |
| 11 completes | Control job 10749750 `COMPLETED` in 2 h 18 min; 19 of 19 stage jobs `COMPLETED`, no `TIMEOUT` or `OUT_OF_MEMORY`; `checkRun.sh` OK; `FRAGLEN=175`. Rates and memory set from `sacct` (tables above) |
| 12 upstream | Identical for all 4 samples: `idxstats` read counts and `@SQ` lines of the postProcess BAMs, uncompressed `clean.bed.gz` (md5), bigWigs (byte for byte), predictd lengths (167, 170, 180, 183) and model files, MACS `broadPeak`, `gappedPeak` and `xls` after renaming `<S>___MD` to `<S>`. `FRAGLEN` 175 equals the baseline's `--extsize`; `macsPairs.txt` and `macsPeaksMerged.saf` (121,100 peaks) identical |
| 13 counts | featureCounts 2.0.6 `-O -Q 10 -p --countReadPairs` on the baseline BAMs and SAF equals the run's 2.1.1 matrix in all 4 x 121,100 cells and in the `.summary`; DESeq2 scale factors recomputed from it equal the run's to 15 digits. Checks A (ratio 1.0000) and B (diff 0) pass for all 4. The v0.9.0 matrix is about twice this one (reads, not fragments): these are the corrected counts for Proj_15235_B |
| 14 postPipe | Exit 0 with the project's `manifest.txt`; QC xlsx and a 4-page PDF (fragments in peaks, fraction, two PCA plots) in `out/qc/`; `___sigPeaks` and `___Volcano` in `out/qc/peaks/` |

Not done (spec section 6): item 17, b37 (`00.ISSUES.md` #1).
