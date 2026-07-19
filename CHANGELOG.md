# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
