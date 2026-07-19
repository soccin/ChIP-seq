# Results for ChIPSeq Pipeline (v0.9.0)

Output folders:

_metrics_:
PICARD metrics from mapping phase

_alignments_:
RAW unprocessed BAM files. Note before using the BAM's we post process them according to the ENCODE recipies which does the following:

- remove unmapped, mate unmapped, non-proper paired reads
- remove reads with MAPQ < 30
- remove duplicates
- remove Failed QC reads

_chipSeq/macs_: output of running `macs2 callpeak` on BED files derived from the post processed bams. The peak type is chosen once for the whole run, not per sample. Broad peaks are the default, for IP's with non-focal/broad binding:

- `--broad`

for IP's with focal transcription factors the pipeline is run with `-n`, which calls the narrow peak set instead (no `--broad`).

In both cases the following arguments are used:

- `-p 0.01` - peak cutoff on the p-value rather than the q-value
- `--nomodel --shift 0 --extsize <FRAGLEN>` - skip MACS model building and extend reads by a fixed fragment length. `<FRAGLEN>` is the median of the per-sample fragment lengths estimated by `macs2 predictd`, so a single value is used for every sample in the run.

See the MACS2 website [https://github.com/taoliu/MACS] for more information on the output.

There is also an aggregate peak file there which combines the peaks from all samples
- macs/macsPeaksMerged.saf
and a file that counts the coverage for each sample within these combined peaks
- macs/peaks_raw_fcCounts.txt

_chipSeq/qc_: QC report/plots for ChIP related QC.

- qcChIPSeq_PROJECT-NUM_.xlsx - is a report for the total number of peaks found and the number of significant peaks. Low number of peaks could indicate an issue with the ChIP.

- qcChIPSeq_PROJECT-NUM_.pdf - plots of the number/percentage of mapped reads that fall in MACS peaks and a PCA plot of the aggregate peak counts for the samples.

_chipSeq/bw_:
Normalized (to 10million reads) bigwig files for loading IP profiles into IGV [http://software.broadinstitute.org/software/igv/].

_chipSeq/annote_:
Annotation of MACS2 peak file using HOMER `annotatePeaks.pl`. See the HOMER website [http://homer.ucsd.edu/homer/ngs/annotation.html] for more details on the output.

_chipSeq/diff_:
Differential analysis of the aggregate peak counts using `edgeR`, run for each requested pairwise group comparison. Only present if a differential analysis was requested.

- `*_DiffPeaksEdgeR_V3.xlsx` - differential peaks per comparison, with significant peaks annotated using ChIPseeker. Significance cutoff is FDR < 0.05.
- `*_DiffPeaks_V3.pdf` - PCA plot plus MA and volcano plots for each comparison.

Comparisons are specified as group pairs; the sign convention is Group2 - Group1, so a positive logFC means higher in Group2.

