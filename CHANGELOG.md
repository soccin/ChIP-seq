# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

To be released as 1.0.0 once the full-size validation in
`00.SPEC_v1.0.0.md` passes. Breaking: Slurm only, counts change, names change.

### Breaking
- IRIS/Slurm only. Submit with `sbatch ChIP-seq/pipe.sh`; LSF support
  removed. Logs in `SLURM.*/`. One run per analysis directory.
  `-o|--outdir` removed; output is always `out/`.
- The count matrix now holds fragment counts. featureCounts 2.0.6 was run
  with `-p` and without `--countReadPairs` from v0.7.0 (2024-05-04) to
  v0.9.0, so every count column was reads (2x). Replaced with subread
  2.1.1 (static binary, md5 `0eff90d8e30311160dbddda542f17383`) and
  `-p --countReadPairs`; single-end runs count reads without `-p`. Counts
  from any v0.7.x to v0.9.x run are not comparable with these.
- Genomes are named by build tag only (`b37`, `hg38`, `mm10`), from
  `lib/genomes/GENOMES.tsv`; GRCh38 (GDC `GRCh38.d1.vd1`) is `hg38`.
  `annotateWithHomer.sh` and `R/diffAnalysisPairwise.R` take no genome
  argument and read it from `00.RUNSTATUS.txt` (`human`/`mouse` removed).
  `mm10+b37` is no longer supported.
- Samples are named by their `@RG SM` tag; outputs are
  `out/<SM>_postProcess.*`; pairing-file and manifest names must equal a
  sample exactly. Unchanged for every past project checked (Proj_16852,
  Proj_15235, Proj_15235_B).
- `lib/genomes/<TAG>.genome` lists autosomes, X and Y only (b37 24, hg38
  24, mm10 21 contigs); it is both the POST region set and the chromosome
  sizes. BigWigs no longer list `chrM`/`MT`, which had no reads.
- `postPipe.sh` takes no arguments and refuses a run that
  `bin/checkRun.sh` does not pass.

### Not tested
- b37 is supported but was not tested on IRIS; no b37 ChIP data was
  available. Use with care until `00.ISSUES.md` #1 is closed.

### Added
- Count checks, from ATAC-seq v2.0.0 (devs 28978c9): the run fails if any
  sample's fragment total is not about half its mapped reads (Check A,
  0.95 to 1.05; reads for `-s`), or if its Assigned total differs from an
  independent samtools count (Check B). Report in
  `out/macs/peaks_raw_fcCounts.checks.txt`. The Count stage is
  `bin/countPeaks.sh`; `bin/checkCounts.sh` also audits old runs whose
  BAMs `postPipe.sh` moved to `out/bam/`.
- Run status: `00.RUNSTATUS.txt`, `SLURM.CTRL/jobs.tsv`, `#CHIP_EXIT=`
  trailers, `bin/checkRun.sh`; `COMPLETED` requires every deliverable;
  walltime classes; `/localscratch` temp files.
- Up-front checks before anything is submitted: the genome of every BAM,
  every BAM header against the genome file, one unique `SM` per BAM,
  paired or single-end reads against `-s`, the pairing file readable.
- `out/pipeline_info/samples.tsv` (SAMPLE, BAM, PBAM); `cmd.sh.log` also
  records fragment length, samples, pairing file, control job and tool
  versions, and is written on the no-pairing-file path too.
- GRCh38 support (`hg38`): genome file, HOMER `hg38`,
  `TxDb.Hsapiens.UCSC.hg38.knownGene`.
- `bin/genomeInfo.sh`, `bin/getGenomeBuildBAM.sh` (from ATAC-seq; adds the
  IRIS `hg38+mm39` hybrid as a named, unsupported build),
  `bin/picardJvm.sh` (from PEMapper), `00.SETUP.sh`, `00.ISSUES.md`,
  `docs/SLURM_PORT.md`.

### Changed
- `bin/picardV2`: Java 23 module (the old lookup loaded Java 1.7 on IRIS),
  heap from the job's memory, `TMP_DIR` the job's node-local `TMPDIR`;
  picard stays 2.25.5.
- POST: intermediates on `/scratch/core001/bic/$USER/ChIPSeq` and removed
  on exit; regions from the genome file (same contigs as the old `@SQ`
  regex for b37, hg38 and mm10); relative `.bam.bai` link in all three
  variants; arguments `GENOME SAMPLE BAM`.
- A `macs2 predictd` failure fails the BW job; the median fragment length
  must be a positive integer.
- A pairing-file control that matches no sample is fatal (it silently
  became "no control" before); `na`, in any case, means no control.
- QC plot titles say fragments: "Fragments in MACS Peaks", "Fraction of
  Fragments in MACS Peaks".
- R scripts find `tools.R` from their own directory and do not depend on
  `~/.Rprofile`; they run with `Rscript --no-init-file`.

### Fixed
- Every error exits non-zero: usage paths, `quit()` in
  `diffAnalysisPairwise.R` and `qc_ChIPSeq_02.R`, `exit -1`.
- `set -e` and `pipefail` in every stage script, with a FATAL line naming
  the failed step (shown in `MESSAGE`).
- `makeBigWigFromBEDZ.sh` read the BED path as the scale factor when given
  three arguments.
- `annotateWithHomer.sh`: genome from the run instead of hardcoded hg19;
  `chr` added only for b37 (mm10 peaks became `chrchr1` before); peak
  contigs checked against the genome; temp file removed; shebang fixed.

### Removed
- To `attic/`: `lsf.sh`, `bin/bSync`, `bin/detectOS.sh`,
  `bedGraphToBigWig`, `bin/bedToBigBed`, `00.SETUP.cmds`, `mergeSamples.sh`,
  the root `getGenomeBuildBAM.sh`, `lib/genomes/mm10+b37.genome`,
  `mouse_mm10.genome`, `mouse_mm10.genome.bed`. `01.TODO.md` is replaced by
  `00.ISSUES.md` and `attic/zzISSUES.md`.
- `merge_pngs()` in `diffAnalysisPairwise.R` (never called; needed
  ImageMagick, which IRIS does not have).

### Not ported
- `deliverResults.sh` (`00.ISSUES.md` #2, `POST-IRIS-PORT.md`).

## [0.9.0] - 2026-07-19

### Added
- Warning when pipe.sh is run without --pairing-file, stating up front that
  MACS peak calling will not run and the pipeline stops after Stage2
- CLAUDE.md with pipeline architecture and working-directory conventions

### Removed
- docs/RESULTS_wDiff.md and its rendered PDF, consolidated into docs/output.md.
  The two were near-duplicate results descriptions that had drifted apart

### Changed
- docs/output.md is now the single results description, following the nf-core
  convention, and gained the differential analysis section from RESULTS_wDiff.md
- Rewrote pipe.sh usage text: per-option descriptions with defaults, pairing
  file format (including the 'na' target-only form), and an example command
- Stage3 message when no pairing file is given now states that MACS is
  skipped, replacing the misleading "not yet implemented" wording

### Fixed
- Corrected MACS arguments in docs/output.md to match callPeaks_ChIPseq.sh:
  cutoff is `-p 0.01` not `-q 0.01`, and `--broad-cutoff` is never passed.
  Documented `--nomodel --shift 0 --extsize` and clarified that peak type is
  chosen once per run rather than per sample
- Corrected post-processing MAPQ filter in docs/output.md from 10 to 30
- Updated stale setup references from CMD.INSTALL.MACS/CMDS.INSTALL.MACS to
  00.SETUP.cmds in pipe.sh and docs/install.md
- Delivery email template in deliverResults.sh referenced a nonexistent
  docs/RESULTS.md under a misspelled directory; now points to ChIP-seq/docs/output.md
- Unified version references across README.md, VERSION.md, docs/output.md,
  and getPackageVersions.R, which had drifted to 0.8.0 and 0.7.x

## [0.8.0] - 2025-11-01

### Added
- Pairwise differential analysis implementation (R/diffAnalysisPairwise.R)
  - Handles pairwise comparisons for ChIP-seq differential peak analysis using edgeR
  - Generates MA plots and volcano plots for each comparison
  - Annotates significant peaks using ChIPseeker
  - Outputs Excel workbook with summary statistics and detailed results
  - Supports both human (hg19) and mouse (mm10) genomes
- Delivery email message template output functionality
- MACS2 installation references in documentation
- Comprehensive usage help with examples

### Changed
- **BREAKING**: Set stricter FDR cutoff from 1.1 (no filtering) to 0.05 for stringent significance threshold
- **BREAKING**: Output filename changed from `DiffPeaksEdgeRv2_NoFilt.xlsx` to `DiffPeaksEdgeRv2.xlsx`
- Update MACS installation to use Python 3.10.2
- Update bedtools module loading to use generic version
- Template formatting with additional newlines for readability
- Moved previous differential analysis scripts to attic/

### Fixed
- Error handling improvements in pipeline
- MACS logging in makeBigWigFromBEDZ.sh - stderr now sent to both log file and stdout for LSF visibility

### Refactored
- Apply clarity-first tidying to differential analysis (R code quality)
- Apply tidyverse style guide to differential analysis (R code formatting)

## [0.7.2] - 2024-05-04

### Added
- v0.7.2 documentation in HTML/PDF format

### Changed
- Freeze for version 0.7.2 release

## [0.7.1] - 2024-05-04

### Added
- Delivery email template

### Changed
- Moved results to docs/output directory for better organization

## [0.7.0] - 2024-05-04

### Changed
- Update featureCounts to latest version
- Move featureCounts executable to bin directory

### Documentation
- Add featureCounts installation information in docs

[Unreleased]: https://github.com/soccin/ChIP-seq/compare/v0.9.0...HEAD
[0.9.0]: https://github.com/soccin/ChIP-seq/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/soccin/ChIP-seq/compare/v0.7.2...v0.8.0
[0.7.2]: https://github.com/soccin/ChIP-seq/compare/v0.7.1...v0.7.2
[0.7.1]: https://github.com/soccin/ChIP-seq/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/soccin/ChIP-seq/releases/tag/v0.7.0
