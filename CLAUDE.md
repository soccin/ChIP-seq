# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Execution environment

IRIS (MSKCC, Slurm) only. The JUNO/LSF path was removed in v1.0.0; nothing
here runs on a laptop. On IRIS the pipeline can be run and tested directly:
`bash -n` and `Rscript --no-init-file -e 'parse(file=...)'` for syntax, the
helper scripts by hand, and test runs with `sbatch` in a fresh directory
(the v1.0.0 test directories are `PortATACseq/PortChIPseq/test*`, inputs in
`PortChIPseq/inputs/`, test scripts in `PortChIPseq/scripts/`).

There is no test suite, linter or CI in this repo (`00.ISSUES.md` #12).

## External dependencies

| What | Where | Version |
| --- | --- | --- |
| samtools | `module load samtools/1.20` (`bin/loadTools.sh`) | 1.20 |
| bedtools | `~/bin/bedtools` on PATH (no module exists) | 2.31.1 |
| Java | `module load java/23.0.1` (`bin/picardV2`) | 23.0.1 |
| Slurm wrappers | `~/bin/sbatch` (exports `SBATCH_SCRIPT_DIR`), `~/bin/bsub` (openlava shim), `~/bin/bSync.sh` | |
| R | on PATH, packages in `getPackageVersions.R` | 4.5.1 |
| HOMER | `/usersoftware/core001/common/RHEL_8/homer`, genomes hg19, hg38, mm10 (`CHIP_HOMER_DIR`) | |
| MACS2 | `venv/`, built by `. 00.SETUP.sh` | 2.2.9.1, numpy 1.23.0 |

Vendored in `bin/`: featureCounts (subread 2.1.1, static, md5
`0eff90d8e30311160dbddda542f17383`), picard 2.25.5
(`bin/jar/2.25.5/picard.jar`, run only through `bin/picardV2`),
`wigToBigWig`.

Nothing depends on `~/.Rprofile`: the R scripts use `length()`, `str_c()`
and so on, and `pipe.sh`/`postPipe.sh` run them with
`Rscript --no-init-file`. Do not reintroduce `len()`, `cc()`, `DATE()` or
`getSDIR()`.

## Working-directory contract

The pipeline is run from an analysis directory that holds the clone (any
name) and gets `out/`:

```
Proj_15235_B/ChIP/          <- cwd for all commands; one run per directory
├── ChIP-seq/               <- this repo
├── 00.RUNSTATUS.txt        <- run record (pipe.sh)
├── SLURM.CTRL/ SLURM.0*/   <- control and stage logs
└── out/                    <- all output; fixed, there is no -o
```

The R scripts find `tools.R` from their own location (`commandArgs()`
`--file=`), so the clone can have any name. `out/` is fixed: `-o` was
removed because the QC and differential scripts and the run records all
assume it (`00.ISSUES.md` #4).

## Running

On hold since 2026-10-12: `pipe.sh` exits 1 at once with a "do NOT use
this pipeline" message (block marked `ON HOLD` after the `#SBATCH`
lines). To test or resume, delete that block on a branch first; see
`00.SPEC_v1.0.0.md` section 0.

```bash
sbatch ./ChIP-seq/pipe.sh --pairing-file PAIRS BAM1 [BAM2 ...]
./ChIP-seq/bin/checkRun.sh          # 0 ok, 1 failed, 2 running, 3 unknown
./ChIP-seq/postPipe.sh              # refuses a run checkRun.sh does not pass
find out/macs | egrep "_peaks.(broad|narrow)Peak$" | xargs -n 1 bsub ... ./ChIP-seq/annotateWithHomer.sh
Rscript --no-init-file ChIP-seq/R/diffAnalysisPairwise.R manifest.csv comps.csv [RUNTAG]
```

`postPipe.sh` moves BAMs to `out/bam/` and BEDs to `out/bed/`, runs QC
stage 1, and QC stage 2 if `manifest.txt` (`SampleID<TAB>Group`, no
header) exists. `annotateWithHomer.sh` and `diffAnalysisPairwise.R` read
the genome from `00.RUNSTATUS.txt`; they take no genome argument.
`deliverResults.sh` is not ported (`POST-IRIS-PORT.md`).

## Architecture

`pipe.sh` is the control job, a copy of the ATAC-seq model (`devs`
28978c9; `docs/SLURM_PORT.md`). It runs under `sbatch` on `cmobic_cpu`
with `--qos=priority`, submits every stage job with `chipSub` (through
`~/bin/bsub` and `bin/runStage.sh`) and blocks between stages with
`waitStage` = `bSync` + `bCheck` (`bin/slurmTools.sh`), which reads job
states from `sacct`. Stage jobs are `SHORT` (1:59, short partitions) or
`LONG` (12 h, qos priority), chosen per job by `runClass` from input size
and the `RATE_*` minutes per GB.

Up front, before any submission, `pipe.sh` checks: the genome of every BAM
(`bin/getGenomeBuildBAM.sh`, all equal, a row of `GENOMES.tsv`), every BAM
header against `lib/genomes/<GENOME>.genome`, one unique `@RG SM` per BAM,
paired or single-end reads against `-s`, and that the pairing file is
readable. It writes `out/pipeline_info/samples.tsv` (SAMPLE, BAM, PBAM),
which every later step reads.

Stages (job name `qChIP-Seq_<STAGE>_<pid>`, log `SLURM.<nn>.<STAGE>/`):

1. **POST** — one of `postMapBamProcessing_ChIPSeq.sh` (default),
   `_NoPP.sh` (`--proper-pair-off`), `_SE.sh` (`-s`), arguments `GENOME
   SAMPLE BAM`. ENCODE3 filter (MAPQ >= 30, `-F 3852`, contigs of the genome
   file), picard sorts, writes `out/<SM>_postProcess.bam` and
   `.clean.bed.gz`. Intermediates under `/scratch/core001/bic/$USER/ChIPSeq`,
   removed on exit.
2. **BW** — `makeBigWigFromBEDZ.sh`: bigWig scaled to 10M reads, and
   `macs2 predictd`, whose `.log` is the only source of fragment length. A
   predictd failure fails the job.
   Then, in the control job, the median fragment length
   (`getMedianFragmentLengthFromPredictDFile.R`), which must be a positive
   integer; one value for every MACS call. Without a pairing file the run
   stops here, `COMPLETED`.
3. **CALLP** — `generateMACSArgs.R` (control job) resolves the pairing file
   against `samples.tsv` by exact name into `out/macsPairs.txt`;
   `callPeaks_ChIPseq.sh` runs one MACS call per line. Broad by default.
4. **MergePeaks** — `mergePeaksToSAF.sh`, peaks within 500 bp, one SAF.
5. **Count** — `bin/countPeaks.sh`: featureCounts 2.1.1
   `-p --countReadPairs` (no `-p` with `-s`), then Check A. Beside it, one
   **CountCheck** job per sample (`bin/assignedCheck.sh`, Check B); the
   control job compares them (`bin/checkCounts.sh`) and writes
   `out/macs/peaks_raw_fcCounts.checks.txt`.
6. **DESEQ** — `getDESeqScaleFactors.R`, the *inverse* of DESeq2's
   `sizeFactors` (deliberate; `NOTES.md`).

Then a `sacct` sweep of every stage job and the deliverables check;
`STATUS=COMPLETED` is written last.

### Run records

- `00.RUNSTATUS.txt`: `STATUS` (RUNNING, COMPLETED, FAILED), `STAGE`,
  `RC`, `MESSAGE` (the reason, from the failing check or the failed job's
  first FATAL log line), `GENOME`, `SAMPLES`, `FRAGLEN`, options, control
  job id. Written by `pipe.sh` and its `EXIT` trap.
- `SLURM.CTRL/jobs.tsv`: every stage job, its sample and log.
- `#CHIP_EXIT=<rc>`: last line of every stage log (`bin/runStage.sh`).
- `out/pipeline_info/cmd.sh.log`: version (`git describe --dirty`), options,
  genome, samples, fragment length, tool versions.

One run per directory: `pipe.sh` refuses a directory with
`00.RUNSTATUS.txt`, `SLURM.CTRL/jobs.tsv`, `out/` or `00_PIPE_DONE`. No
resume.

### Temp files

Never `/tmp`. `bin/runStage.sh` gives each stage job a `TMPDIR` under
`/localscratch/$USER` (`CHIP_LOCAL_TMP`) and removes it; picard, sort,
MACS2, featureCounts and R use it.

## Genome builds

A genome is named by its build tag everywhere: `b37` (GRCh37, contigs
`1`..`22`, `X`, `Y`, `MT`), `hg38` (GRCh38, `chr` contigs; the GDC
`GRCh38.d1.vd1` reference; ATAC-seq still calls it `b38`) and `mm10`.
Species words and stand-in names are errors.

- `lib/genomes/GENOMES.tsv` holds every genome-dependent value (MACS
  genome size, HOMER genome, TxDb, OrgDb, whether to add `chr` for
  annotation). A supported genome is a row of it. Shell reads it with
  `bin/genomeInfo.sh TAG FIELD`, R with `genomeInfo()` in `tools.R`; both
  fail on an unknown tag or field. No stage script names a genome.
- `lib/genomes/<TAG>.genome` lists autosomes, X and Y only. It is the POST
  region set and the chromosome sizes. For `samtools view -L` write the
  3-column BED (`contig 0 length`) first, as the POST scripts do: samtools
  reads the 2-column file as `contig position` and selects one base per
  line. Mitochondrial reads are dropped before peak calling, unlike
  ENCODE's ChIP pipeline (`00.ISSUES.md` #14).
- Nobody types the genome of a run: `pipe.sh` detects it and records it;
  post-run tools read `GENOME=` from `00.RUNSTATUS.txt`.
- `bin/getGenomeBuildBAM.sh` knows more builds than are supported, so that
  unsupported ones get a named error. ATAC-seq copied it (branch
  `fix/genome-tags`); keep the two byte-identical.
- `bin/genomeInfo.sh` and `GENOMES.tsv` differ from ATAC-seq's: here the
  table has `HOMER` and the script has `--tags` and the `GENOME_FILE`
  field; ATAC's table adds `MITO`. Shared columns have the same names and
  values.
- b37 is supported but untested on IRIS (`00.ISSUES.md` #1).

Adding a build: a row in `GENOMES.tsv`, `lib/genomes/<TAG>.genome`
(checked against the reference `.fai`), the md5 entry if missing, the HOMER
genome and annotation packages, and a test run.

## Reliability rules

The featureCounts bug (v0.7.0 to v0.9.0 counted reads, not fragments, in
every column, and nothing noticed for two and a half years) is why these
exist:

1. Errors exit non-zero with a message that names the cause (`exit 1`,
   `quit(status = 1)`, `fatal()`); usage errors too.
2. No silent fallback: a missing input, an unmatched name, a missing
   fragment length, an unknown genome or a failed check stops the run.
3. Names are matched exactly: samples by `SM`, pairing-file and manifest
   entries by string equality, genomes by tag. No regex or prefix matching
   on names.
4. The run detects what it can and checks what it is told (`-s` against
   the reads, BAM headers against the genome file).
5. Counts are never accepted unchecked (Checks A and B); `COMPLETED`
   requires every deliverable; `postPipe.sh` needs `checkRun.sh` to pass.
   No Bioconductor package (Rsubread, `GenomicAlignments`, ...) in a
   pipeline step that produces counts: their APIs change without warning.
   The standalone subread binary is not Bioconductor.
6. No dependence on personal configuration.
7. Every run records what it used.

## Conventions

- Shared with ATAC-seq: `bin/slurmTools.sh`, `bin/runStage.sh`,
  `bin/loadTools.sh`, `bin/checkRun.sh`, `bin/countPeaks.sh`,
  `bin/checkCounts.sh`, `bin/assignedCheck.sh` were copied from ATAC-seq
  `devs` 28978c9 with `ATAC` renamed `CHIP`. Local changes: `loadTools.sh`
  is safe under `set -u`; `countPeaks.sh` and `checkCounts.sh` take
  `--single-end`, and `countPeaks.sh` puts the first failing sample into
  its FATAL line (so `MESSAGE` names it); `checkCounts.sh` finds BAMs
  moved to `out/bam/`. Diff again when ATAC releases v2.0.0 and offer
  these back. `bin/getGenomeBuildBAM.sh` went the other way: ATAC-seq
  took this repo's copy.
- `bin/picardJvm.sh` is PEMapper's (v_4.1.0) with `PEMAP_` renamed `CHIP_`.
- Version lives in five places that must move together: `README.md`,
  `VERSION.md`, `CHANGELOG.md` (entry and compare links), `docs/output.md`
  and the `cat()` banner in `getPackageVersions.R`. `pipe.sh` reports
  `git describe`, so the tag must be pushed.
- `docs/output.html` and `.pdf` are rendered from `docs/output.md` with
  `make -C ChIP-seq/docs` (pandoc, xelatex) and must be regenerated when it
  changes.
- Superseded scripts go to `attic/` rather than being deleted. Open issues
  in `00.ISSUES.md`, closed ones in `attic/zzISSUES.md`.
- Shell: `set -e -o pipefail` (and `-u` where safe), `SDIR` self-location,
  `|| { echo "FATAL ERROR: ..."; exit 1; }` on every step that can fail;
  `logReason` puts the first FATAL line of a failed job into `MESSAGE`.
