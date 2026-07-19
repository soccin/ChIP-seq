# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Execution environment

This repo is edited on macOS but only runs on the MSKCC Juno Linux cluster. `module`,
`bsub`/LSF, `picardV2`, and the MACS2 `venv` do not exist locally, so `pipe.sh` and the
stage scripts cannot be executed here — a local run dies at `module load bedtools` under
`set -eu` and tells you nothing.

Verify shell changes with `bash -n <script>`, or by rendering a function/heredoc fragment
in isolation. Note that a fragment extracted with `sed`/`head` may report a spurious
"unexpected end of file" from the truncation itself; trust `bash -n` on the whole file.
Hand anything needing a real run back to the user.

There is no test suite, linter, or CI in this repo.

## Working-directory contract

The pipeline is not run from inside this repo. It is cloned as a subdirectory named
`ChIP-seq` inside a per-project analysis directory, and every command is issued from that
parent:

```
Proj_15235_B/ChIP/          <- cwd for all commands
├── ChIP-seq/               <- this repo (name matters)
├── pipeline/               <- symlink to delivery dir with alignments/
└── out/                    <- ODIR
```

This name is load-bearing: `qc_ChIPSeq_01.R` and `qc_ChIPSeq_02.R` do
`source("ChIP-seq/tools.R")` with a hardcoded relative path.

Equally load-bearing is `out`. Although `pipe.sh` takes `-o|--outdir`, the QC and
differential scripts hardcode `out/macs` (`qc_ChIPSeq_01.R:9`, `qc_ChIPSeq_02.R:28,67`,
`R/diffAnalysisPairwise.R:177`). Changing `-o` silently breaks every post-pipeline stage.

The R scripts also call `len()`, which is not defined anywhere in this repo — it comes
from the user's `~/.Rprofile`. Read that file before assuming an undefined R helper is a bug.

## Running

Setup (once per clone) — must be sourced, not run, and is the only way `venv/` appears:

```bash
. 00.SETUP.cmds
```

Full run, from the project dir (see `CMDS.RunProjectTMPL`):

```bash
bsub -o LSF.00.CTRL/ -W 359 ./ChIP-seq/pipe.sh \
    --pairing-file pipeline/*_sample_pairing.txt \
    pipeline/alignments/*bam

./ChIP-seq/postPipe.sh out
find out/macs | egrep "broad|narrow" | xargs -n 1 ./ChIP-seq/annotateWithHomer.sh
./ChIP-seq/deliverResults.sh <RESDIR>
```

`postPipe.sh` runs QC stage 1, then runs QC stage 2 only if `manifest.txt`
(`SampleID<TAB>Group`, no header) exists in the cwd. Otherwise it prints the command to
run by hand.

Differential analysis is separate and manual:

```bash
Rscript --no-save ChIP-seq/R/diffAnalysisPairwise.R <human|mouse> manifest.csv comps.csv [RUNTAG]
```

## Architecture

`pipe.sh` is a synchronous LSF driver, not a workflow engine. Each stage fans work out
with `xargs -n N bsub ...` under a job tag containing `$$`, then blocks on
`bSync <TAG>` (`bin/bSync`, which polls `bjobs` and submits a `-K` dependent job) before
the next stage starts. `lsf.sh` wraps `bsub` purely to echo each submission. Stages 4-6
instead chain with `-w "post_done(...)"` and only the final DESeq job is `bSync`'d.

Consequences worth knowing: the driver must stay alive for the whole run, so it is itself
submitted with a long wall clock; and stage boundaries are the natural place to resume a
failed run (there are commented-out `if [ "" ]; then` guards in `pipe.sh` used to skip
stages 1-2 during debugging).

Stages:

1. **Post-processing** — one of three sibling scripts by read-type flags:
   `postMapBamProcessing_ChIPSeq.sh` (default), `_NoPP.sh` (`--proper-pair-off`, for
   samples where translocations matter), `_SE.sh` (`-s`). Applies the ENCODE3 filter
   recipe (MAPQ>=30, `-F 3852`, main chromosomes only, no MT/decoys), then emits both a
   filtered `*_postProcess.bam` and a `*.clean.bed.gz`.
2. **bigWig** — `makeBigWigFromBEDZ.sh`, scaled to 10M reads unless a scale factor is
   passed. Also runs `macs2 predictd` per sample; the resulting `.log` files are the
   only source of fragment length.
3. **Peak calling** — `generateMACSArgs.R` resolves the pairing file against actual BAM
   names into target/control pairs, which `callPeaks_ChIPseq.sh` consumes two at a time.
   Broad peaks are the default; `-n` switches to narrow.
4. **Merge peaks** — `mergePeaksToSAF.sh`, merging within 500bp into a single SAF.
5. **Counts** — `featureCounts` over the merged SAF against all post-processed BAMs.
6. **Scale factors** — `getDESeqScaleFactors.R`, writing the *inverse* of DESeq2's
   `sizeFactors` (this inversion is deliberate; see `NOTES.md`).

The fragment length passed to MACS is a single scalar: the median across all samples of
the per-sample `predictd` estimates, computed by
`getMedianFragmentLengthFromPredictDFile.R`. One outlier sample shifts every peak call.

Provenance is written to `$ODIR/pipeline_info/cmd.sh.log` (git tag, remote, cwd, args) at
the end of every run, and `date` is written to `00_PIPE_DONE` / `01_POST_DONE` as stage
sentinels.

### Pairing file

Two columns, control then target. `na` in column 1 means target-only, which
`generateMACSArgs.R` maps to the sentinel `_NoCTRL` and `callPeaks_ChIPseq.sh` turns into
an omitted `-c`. Without `--pairing-file` the pipeline warns up front and stops after
stage 2.

Sample matching is by regex (`<name>_.*postProcess`), not exact string. A name that is a
prefix of another sample's name matches both and aborts as "degenerate target/control".

### Genome handling

`getGenomeBuildBAM.sh` identifies the build by md5 of the sorted `@SQ` header lines, and
recognizes far more builds than the pipeline actually supports. Downstream,
`callPeaks_ChIPseq.sh` and `makeBigWigFromBEDZ.sh` accept only `mm10*` and `b37`, and
`lib/genomes/` holds only b37, mm10, mm10+b37, mouse_mm10. A BAM that resolves to `hg19`,
`GRCh37-lite`, or `mm9Full` passes the initial check and then fails several stages in.
Adding a build means touching the md5 case, both downstream case statements, and
`lib/genomes/`.

`annotateWithHomer.sh` is hardcoded to hg19 regardless of the detected genome, and prints
a warning saying so. It also depends on a HOMER install at a hardcoded path.

## Conventions

- Version lives in five places that must move together, and has drifted before:
  `README.md`, `VERSION.md` (which also pins the full R/Python package set),
  `CHANGELOG.md` (entry plus the compare links at the bottom), `docs/output.md`, and the
  `cat()` banner in `getPackageVersions.R`. `pipe.sh` reports its own version from
  `git describe`, so it needs no edit but does need the release tag to be pushed.
- `docs/output.html` and `docs/output.pdf` are rendered from `docs/output.md` and must be
  regenerated when it changes; the render command is not captured in the repo.
- Superseded scripts go to `attic/` rather than being deleted.
- Shell: `set -eu`, `SDIR="$( cd "$( dirname "$0" )" && pwd )"` for self-location, and
  `bin/` prepended to `PATH` for vendored binaries (`featureCounts`, `wigToBigWig`,
  `picardV2`, `bSync`).
