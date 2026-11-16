# Results for ChIPSeq Pipeline (v1.0.0)

This document describes the files in the `chipSeq` folder of a delivery
and how they were made. The options of your run (genome, peak type,
fragment length, samples) are recorded in
`chipSeq/pipeline_info/cmd.sh.log`.

## Folders

| Folder | Contents | Present |
| --- | --- | --- |
| `bw` | Coverage profiles (bigWig), one per sample | always |
| `macs` | MACS2 peak calls, merged peak set, count matrix | if peaks were called |
| `qc` | QC report, significant peaks and plots per peak call | if peaks were called |
| `annote` | HOMER annotation of the peaks | if requested |
| `diff` | Differential peak analysis (edgeR) | if requested |
| `pipeline_info` | Record of the run | always |
| `docs` | This document and the methods, as PDF and HTML | always |

Sample names are the sample names of the mapped BAM files (the `SM` tag
of the read group).

The delivery can also have the folders `bams` (the mapped BAM files) and
`metrics` (alignment metrics) beside `chipSeq`. They come from the
mapping pipeline and are not described here.

## Processing

### Read filtering

The mapped BAM files are filtered before any other step, following the
ENCODE ChIP-seq pipeline. A read is kept only if:

- it is mapped, with mapping quality (MAPQ) of at least 30
- it is a primary alignment (no secondary or supplementary alignments)
- it is not marked as a duplicate and did not fail the sequencer QC
- it is on an autosome, chrX or chrY; reads on the mitochondrial genome
  and on unplaced, random, decoy and viral contigs are removed
- for paired-end data, its mate is mapped and passes the same filters,
  and the two reads are a proper pair (same chromosome, expected
  orientation and distance)

Two options change the last rule; `cmd.sh.log` shows which were used:

- `PROPER_PAIR: No` - the proper pair requirement was dropped, to keep
  pairs that span a rearrangement. Both reads still have to be mapped.
- `SE: Yes` - single-end data; there are no mates and no pair rules.

The filtered BAM files are not delivered. The BAM files in `bams`, when
present, are the mapped files before this filtering.

### Fragment length

`macs2 predictd` estimates the fragment length of each sample from its
filtered reads. The median over the samples (`FRAGLEN` in `cmd.sh.log`)
is used for every peak call of the run.

## _chipSeq/bw_

`<SAMPLE>_postProcess.clean.10mNorm.bw`: coverage of the filtered reads
of the sample, scaled to 10 million reads, so that samples of different
sequencing depth can be compared. Both reads of a pair are counted and
reads are not extended. The files can be loaded into IGV
[<https://igv.org>] or any other genome browser.

## _chipSeq/macs_

Peaks are called with `macs2 callpeak` on the filtered reads, once per
line of the pairing file, so once per target sample. Each call has its
own folder, named for the target and its control:

- `<TARGET>_<CONTROL>` - target called against a control (input or IgG)
- `<TARGET>____NoCTRL` - target called without a control

The peak type is chosen once for the whole run (`PEAK_TYPE` in
`cmd.sh.log`):

- `broad` (default, `--broad`) - for histone marks and other factors
  with broad binding. Files: `_peaks.broadPeak`, `_peaks.gappedPeak`,
  `_peaks.xls`.
- `narrow` - for transcription factors and other focal binding. Files:
  `_peaks.narrowPeak`, `_summits.bed`, `_peaks.xls`.

Both use these arguments:

- `-p 0.01` - peaks are called at p-value < 0.01, a permissive cutoff.
  Each peak also has a q-value; the QC step below uses it to select the
  significant peaks.
- `--nomodel --shift 0 --extsize <FRAGLEN>` - MACS2 does not build its
  own model; reads are extended to the fragment length described above.
- `-g` - the effective genome size of the build (`hs` or `mm`).

The `_peaks.xls` file is a tab-delimited text file, not an Excel file.
See the MACS2 documentation [<https://github.com/macs3-project/MACS>]
for the format of the output files.

The other files in the folder are for the run as a whole:

- `macsPeaksMerged.saf` - the merged peak set: the peaks of all calls,
  with peaks closer than 500 bp joined into one. Columns: peak name
  (`Peak_N`), chromosome, start, end, strand. The start is the 0-based
  start of the interval, one base lower than SAF defines it.
- `peaks_raw_fcCounts.txt` - the raw count matrix: for each merged peak
  and each sample, the number of fragments that overlap the peak
  (featureCounts; read pairs for paired-end data, reads for single-end
  data). A fragment that overlaps two peaks is counted in both. The
  count columns are named for the filtered BAM files,
  `out/<SAMPLE>_postProcess.bam`. This is the input for any
  normalization or differential analysis.
- `peaks_raw_fcCounts.txt.summary` - for each sample, the number of
  fragments in peaks (`Assigned`) and outside them.
- `peaks_raw_fcCounts.checks.txt` - two checks of the count matrix
  against independent counts made with samtools. A run is delivered only
  if every line is `OK`.
- `scaleFactorsDESeq2.csv` - one scale factor per sample from the count
  matrix: the inverse of the DESeq2 size factor. Multiply the counts of
  a sample by it to compare samples. The bigWig files do not use it.

In versions v0.7.0 to v0.9.0 `peaks_raw_fcCounts.txt` counted each read
of a pair separately, so every count was twice the number of fragments.
Since v1.0.0 it counts fragments.

## _chipSeq/qc_

- `qcChIPSeq_<PROJECT>.xlsx` - one row per peak call: the number of
  peaks (`TotalPeaks`), the number of significant peaks, those with a
  MACS2 q-value < 0.05 (`SigPeaks`), and their fraction (`PCT`). A low
  number of peaks can indicate a problem with the ChIP.
- `qcChIPSeq_<PROJECT>.pdf`
    - the number of fragments of each sample inside and outside the
      merged peaks
    - the same as a percentage
    - a PCA of the samples from the count matrix (edgeR: low-count peaks
      removed, TMM normalization, log CPM), colored by sample group,
      without and with sample labels. Use it to look for outliers and to
      see whether replicates group together. Samples marked as excluded
      are left out.
- `peaks/` - for each peak call:
    - `<CALL>___sigPeaks__FDR_0.05.xlsx` - the peaks with q-value < 0.05,
      most significant first, with the columns of the MACS2 `_peaks.xls`
      file. Not written when the call has no significant peak.
    - `<CALL>___Volcano__FDR_0.05.png` - fold enrichment against p-value
      for every peak of the call; significant peaks in red.

## _chipSeq/annote_

Annotation of each MACS2 peak file with HOMER `annotatePeaks.pl`: for
each peak the nearest gene, the distance to its TSS and the genomic
feature (promoter, exon, intron, intergenic, ...).

- `<CALL>_peaks.broadPeak_HOMERAnnote.txt` (or `narrowPeak`) - the HOMER
  output
- `<CALL>_peaks.broadPeak_HOMERAnnote.xlsx` - the same table as an Excel
  file

See the HOMER documentation
[<http://homer.ucsd.edu/homer/ngs/annotation.html>] for the columns.

## _chipSeq/diff_

Differential analysis of the count matrix with edgeR, for each requested
comparison of two sample groups. For each comparison only the samples of
the two groups are used: low-count peaks are removed (`filterByExpr`),
the counts are normalized (TMM), and each peak is tested with a
likelihood ratio test (`glmLRT`).

The sign convention is Group2 - Group1. Sheets and plots are named
`comp_group<Group2>-comp_group<Group1>`, the group that is higher at
positive `logFC` first: in `comp_groupKO-comp_groupWT` a positive `logFC`
means the peak is higher in KO.

- `<PROJECT>__DiffPeaksEdgeR_V3.xlsx`
    - `Summary` - the number of significant peaks for each comparison
    - one sheet per comparison with its significant peaks: `logFC`
      (log2), `logCPM`, `PValue`, `FDR`, the position and length of the
      peak, and the nearest gene (`SYMBOL`, `GENENAME`), genomic feature
      (`annotation`) and distance to the TSS from ChIPseeker.
- `<PROJECT>__DiffPeaks_V3.pdf`
    - a PCA of the samples in the analysis, without and with labels
    - for each comparison an MA plot and a volcano plot; significant
      peaks in red

A peak is significant when its `FDR` value is below 0.05. The `FDR`
column holds the p-value adjusted for multiple testing with the Holm
method, which is stricter than the Benjamini-Hochberg false discovery
rate. Only significant peaks are listed in the Excel file.

If there is more than one analysis of a run, the file names have a tag
between the project and `DiffPeaks`.

## _chipSeq/pipeline_info_

The record of the run.

- `cmd.sh.log` - pipeline version, options (proper pairs, peak type,
  single-end), pairing file, genome, fragment length, samples, input BAM
  files, versions of the tools, and the command line
- `samples.tsv` - for each sample, its input BAM file and its filtered
  BAM file

## _chipSeq/docs_

- `output.pdf`, `output.html` - this document
- `methods.pdf`, `methods.html` - the methods of the pipeline as text
  for a paper, with the differences from the ENCODE pipeline, the tool
  versions and the references

## Software

The versions used for a run are in `cmd.sh.log`. As of v1.0.0: samtools
1.20, bedtools 2.31.1, picard 2.25.5, MACS2 2.2.9.1, featureCounts
(subread) 2.1.1, R 4.5.1 with edgeR, DESeq2 and ChIPseeker, and HOMER.

Supported genome builds: `b37` (GRCh37), `hg38` (GRCh38), `mm10`.

Source and change history: <https://github.com/soccin/ChIP-seq>
