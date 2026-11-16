# ChIP-seq Pipeline

### Version 1.0.0

## Overview

Samples were processed with the BIC ChIP-seq pipeline, which takes
aligned, duplicate-marked BAM files and produces filtered alignments,
normalized signal tracks (bigWig), MACS2 peak calls for each target
sample (against its control when there is one), a merged peak set with
a raw fragment count matrix, QC reports and, when requested, an
annotation of the peaks and a differential binding analysis. The code
is at
[https://github.com/soccin/ChIP-seq/releases/tag/v1.0.0](https://github.com/soccin/ChIP-seq/releases/tag/v1.0.0).

Read filtering and peak calling follow the ENCODE 3 ChIP-seq pipeline
(ENCODE-DCC/chip-seq-pipeline, `filter_qc.py` and `macs2.py`) with the
differences described below: where duplicates are marked, the removal
of supplementary alignments, the contigs that are kept, how the
fragment length is estimated, and the duplicate setting of MACS2. The
pipeline does not filter peaks against the ENCODE blacklist and does
not run ENCODE's replicate analysis (IDR, pooled and pseudo-replicate
peaks). The signal tracks, the merged peak set and the count
normalization follow the methods of R. Koche (MSKCC Center for
Epigenetics Research).

Reads were aligned to one of three reference genomes: GRCh37 (b37),
GRCh38 (the GDC `GRCh38.d1.vd1` build, hg38) or GRCm38 (mm10). The
build of a project is recorded in the BAM headers and in
`pipeline_info/cmd.sh.log`, which also holds the options of the run
(peak type, proper pairs, single-end, fragment length). Peaks, signal
tracks and counts use the autosomes, X and Y only.

## Methods

**Alignment.** Reads were aligned with BWA-MEM 0.7.19 (Li 2013) using
the PEMapper pipeline
([https://github.com/soccin/PEMapper](https://github.com/soccin/PEMapper)),
which trims adapters, aligns with `bwa mem -M`, merges lanes, sorts,
and marks duplicates with Picard MarkDuplicates. Duplicates are marked
in the BAM, not removed. The PEMapper and BWA versions are in the `@PG`
header lines of each BAM.

**Read filtering.** With samtools 1.20 (Danecek 2021) each BAM was
reduced to properly paired reads with mapping quality 30 or more,
without unmapped reads, reads whose mate is unmapped, secondary and
supplementary alignments, reads failing vendor QC or duplicates (`-q 30
-f 2 -F 3852`), on the autosomes, X and Y. The reads were then
name-sorted and passed through `samtools fixmate -r`, and the flag
filter was applied again (`-f 2 -F 1804`), as in the ENCODE pipeline,
which removes a read whose mate did not pass the filters above. The
result was coordinate-sorted with Picard SortSam 2.25.5; this filtered
BAM is the input to every later step. Three points differ from ENCODE.
The duplicates removed here were marked by PEMapper on all alignments,
before this filter; ENCODE marks duplicates after its own filter.
Supplementary alignments are removed (flag 2048), which ENCODE's flag
1804 does not include. Reads on the mitochondrial contig and on
unplaced, random, decoy and viral contigs are removed at this step, so
they are not part of the read totals MACS2 uses; ENCODE's ChIP-seq
pipeline keeps them through peak calling. When a run is made without
the proper pair requirement (`PROPER_PAIR: No`, for samples where
rearrangements matter) the `-f 2` condition is dropped and both reads
of a pair still have to be mapped; for single-end data (`SE: Yes`)
there are no mate conditions and no fixmate step.

**Signal tracks.** Each read of the filtered BAM was converted to a BED
interval with bedtools 2.31.1 (Quinlan 2010). For each sample the
coverage of these intervals was computed with bedtools genomecov, each
read as its own interval with no extension, scaled to 10 million reads
(the scale factor is 10 million divided by the number of filtered reads
of the sample), and written as a bigWig with the UCSC wigToBigWig tool
(Kent 2010). Tracks are therefore comparable between samples by
sequencing depth; the DESeq2 scale factors below are not applied to
them.

**Fragment length.** The fragment length of each sample was estimated
from its filtered reads with `macs2 predictd`, and the median over the
samples of the run (rounded down) was used as the one fragment length
of every peak call. ENCODE instead takes the fragment length of each
experiment from a strand cross-correlation analysis. A run stops if
predictd gives no estimate for a sample.

**Peak calling.** Peaks were called with MACS2 2.2.9.1 (Zhang 2008)
from the BED intervals of the reads, once for each target sample,
against its control sample (input or IgG) when the experiment has one
and without a control otherwise. Both reads of a pair are used as
separate reads. The settings are those of the ENCODE pipeline: `-p 0.01
--nomodel --shift 0 --extsize <fragment length>` and the genome size
`hs` or `mm`, with `--broad` for broad peaks (the default; histone
marks and other broad binding) and without it for narrow peaks
(transcription factors and other focal binding). The peak type is one
choice for the whole run. Each read is extended to the fragment length
from its 5' end, and peaks are called at a p-value of 0.01, a
permissive cutoff chosen by ENCODE so that peaks can be ranked
downstream. One setting differs from ENCODE: MACS2's default duplicate
handling was used (`--keep-dup 1`, at most one read per position and
strand) and not `--keep-dup all`. MACS2's output (broadPeak and
gappedPeak, or narrowPeak and summits, and the `_peaks.xls` table) is
delivered unchanged.

**Merged peak set.** The peaks of all calls (broadPeak or narrowPeak
files) were pooled and merged with bedtools merge, joining peaks within
500 bp of each other into one region (R. Koche's method; ENCODE does
not merge peaks). The result is one set of regions shared by every
sample of the project, written in SAF format with the regions named
`Peak_1`, `Peak_2`, ... in genome order. The start of each region is
the 0-based start from bedtools written unchanged into the 1-based SAF,
so each counted region begins one base before the merged peak. Because
peaks are merged and pooled over samples, the number of regions grows
with the number of samples and is not comparable to per-sample peak
counts.

**Counts.** Fragments were counted over the merged peak set with
featureCounts 2.1.1 (Liao 2014) from the filtered BAMs, counting read
pairs as single fragments (`-p --countReadPairs`; reads for single-end
data), with a fragment that overlaps two regions counted in both
(`-O`). Control samples are counted as well. The counts were verified
by two independent checks before use: each sample's fragment total must
be half its mapped reads, within 5%, and the number of fragments in the
merged peaks, counted a second time with samtools alone, must equal
featureCounts' assigned total exactly. A run whose counts fail either
check is stopped. The result of both checks for each sample is
delivered with the counts.

**Normalization.** Size factors were computed from the raw count matrix
with DESeq2 (Love 2014; median-of-ratios method, Anders 2010) and
written as scale factors, the reciprocal of the DESeq2 size factor, so
that a normalized count is the raw count multiplied by the sample's
scale factor. They are delivered with the counts but are not applied to
any delivered file.

**QC.** For each peak call the peaks with a MACS2 q-value below 0.05
were taken as significant; their number and the total number of peaks
are reported, the significant peaks are written to a table, and every
peak is plotted as fold enrichment against p-value. The fraction of
fragments in peaks is the number of fragments featureCounts assigned to
the merged peak set divided by all fragments of the filtered BAM, taken
from the featureCounts summary. It measures signal against background
as ENCODE's FRiP does, but is a different quantity: the regions are the
merged peaks of all samples rather than the sample's own peaks, so it
is not comparable to ENCODE's FRiP cutoffs. A principal component
analysis of the samples was computed from the count matrix with edgeR
(regions with too few counts removed with filterByExpr, TMM
normalization, log2 CPM, not scaled) and plotted colored by sample
group.

**Peak annotation.** When requested, each MACS2 peak file was annotated
with HOMER v5.1 `annotatePeaks.pl` (Heinz 2010) and the HOMER genome of
the build (hg19 for b37, with `chr` added to the contig names; hg38;
mm10), which gives for each peak the nearest TSS, its gene, the
distance and the genomic feature.

**Differential binding.** When requested, pairwise comparisons between
sample groups were tested with edgeR (Robinson 2010; McCarthy 2012).
For each comparison, the samples of its two groups were taken from the
count matrix, regions with too few counts removed with filterByExpr
(over all samples of the analysis, then over the samples of the
comparison), library sizes normalized with TMM, dispersions estimated,
and a negative binomial generalized linear model fitted and tested with
the likelihood ratio test. P-values were adjusted with the Holm method
(Holm 1979; the column is labeled `FDR` in the output) and regions with
an adjusted p-value below 0.05 were called significant. The Holm method
controls the family-wise error rate and is stricter than the
Benjamini-Hochberg false discovery rate. The sign of the log fold
change is group 2 minus group 1 of the comparison: a positive value
means more signal in group 2. Significant regions were annotated with
the nearest gene and genomic feature using ChIPseeker (Yu 2015) with
the UCSC knownGene annotation of the build and a TSS region of 5 kb on
each side. For each comparison, every tested region was plotted as an
MA plot (log fold change against mean log2 CPM) and a volcano plot (log
fold change against p-value, on a log scale), with the significant
regions in red. Samples can be left out of the analysis by group.

## Tool versions

The versions used for a run are recorded in `pipeline_info/cmd.sh.log`.

| Tool | Version |
| --- | --- |
| BWA-MEM (PEMapper) | 0.7.19 |
| samtools | 1.20 |
| Picard SortSam | 2.25.5 |
| bedtools | 2.31.1 |
| UCSC wigToBigWig | v 4 |
| MACS2 | 2.2.9.1 |
| featureCounts (subread) | 2.1.1 |
| HOMER | v5.1 |
| R | 4.5.1 |
| DESeq2 | 1.50.2 |
| edgeR | 4.8.2 |
| ChIPseeker | 1.46.1 |
| TxDb.Hsapiens.UCSC.hg38.knownGene (hg38) | 3.22.0 |
| TxDb.Hsapiens.UCSC.hg19.knownGene (b37) | 3.22.1 |
| TxDb.Mmusculus.UCSC.mm10.knownGene (mm10) | 3.10.0 |
| org.Hs.eg.db, org.Mm.eg.db | 3.22.0 |

## References

Anders S, Huber W. Differential expression analysis for sequence count
data. Genome Biol 2010;11:R106. doi:10.1186/gb-2010-11-10-r106

Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience
2021;10:giab008. doi:10.1093/gigascience/giab008

ENCODE ChIP-seq pipeline. ENCODE Data Coordination Center.
[https://github.com/ENCODE-DCC/chip-seq-pipeline](https://github.com/ENCODE-DCC/chip-seq-pipeline).
Data standards: [https://www.encodeproject.org/chip-seq/](https://www.encodeproject.org/chip-seq/)

Heinz S, Benner C, Spann N, et al. Simple combinations of
lineage-determining transcription factors prime cis-regulatory elements
required for macrophage and B cell identities. Mol Cell
2010;38:576-589. doi:10.1016/j.molcel.2010.05.004

Holm S. A simple sequentially rejective multiple test procedure. Scand
J Stat 1979;6:65-70.

Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and
BigBed: enabling browsing of large distributed datasets. Bioinformatics
2010;26:2204-2207. doi:10.1093/bioinformatics/btq351

Li H. Aligning sequence reads, clone sequences and assembly contigs with
BWA-MEM. arXiv 2013;1303.3997.

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose
program for assigning sequence reads to genomic features.
Bioinformatics 2014;30:923-930. doi:10.1093/bioinformatics/btt656

Love MI, Huber W, Anders S. Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. Genome Biol 2014;15:550.
doi:10.1186/s13059-014-0550-8

McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of
multifactor RNA-Seq experiments with respect to biological variation.
Nucleic Acids Res 2012;40:4288-4297. doi:10.1093/nar/gks042

Picard Toolkit. Broad Institute.
[http://broadinstitute.github.io/picard/](http://broadinstitute.github.io/picard/)

Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for
comparing genomic features. Bioinformatics 2010;26:841-842.
doi:10.1093/bioinformatics/btq033

Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
differential expression analysis of digital gene expression data.
Bioinformatics 2010;26:139-140. doi:10.1093/bioinformatics/btp616

Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP
peak annotation, comparison and visualization. Bioinformatics
2015;31:2382-2383. doi:10.1093/bioinformatics/btv145

Zhang Y, et al. Model-based analysis of ChIP-Seq (MACS). Genome Biol
2008;9:R137. doi:10.1186/gb-2008-9-9-r137
