#
# qc_ChIPSeq_02.R MANIFEST [RUNTAG]
#
# QC stage 2, run by postPipe.sh from the analysis directory with
# manifest.csv. MANIFEST is the CSV that R/diffAnalysisPairwise.R also
# reads: header MapID,SampleID,Group (read_manifest() in tools.R). It
# must list every sample of the run; SampleID is the name in the plots,
# Group colors the PCA, and groups starting with EXC are left out of the
# PCA.
#
# Writes out/qc/qcChIPSeq_<Proj>[_<RUNTAG>].pdf: fragments in and outside
# the merged MACS peaks per sample, and a PCA of the peak counts.
#

args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
library(tidyverse)
library(scales)
library(patchwork)
library(edgeR)
library(ggrepel)
library(pals)
})

SDIR <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |> str_remove("^--file=") |> dirname()
source(file.path(SDIR, "tools.R"))

if(length(args)<1 || length(args)>2) {
    fatal("usage: qc_ChIPSeq_02.R MANIFEST [RUNTAG]")
}
RUNTAG <- if (length(args) == 2) args[2] else ""

samples <- read_samples()
manifest <- read_manifest(args[1], samples, require_all = TRUE)

#
# Sample name (MapID) to the name in the plots (SampleID)
#
sample_id <- set_names(manifest$SampleID, manifest$MapID)

ds=read_tsv("out/macs/peaks_raw_fcCounts.txt.summary",show_col_types=FALSE)
#
# Columns are the postProcess BAM paths; samples.tsv maps them to sample
# names exactly.
#
colnames(ds)=c("Status",pbam_to_sample(colnames(ds)[-1],samples))

ds <- ds |>
    gather(Sample, Count, -Status) |>
    mutate(Sample = unname(sample_id[Sample])) |>
    mutate(Status = str_remove(Status, "_.*")) |>
    group_by(Sample, Status) |>
    summarize(Counts = sum(Count), .groups = "drop") |>
    mutate(Status = ifelse(Status == "Assigned", "InPeaks", "Outside"))

pg0=ggplot(ds,aes(Sample,Counts,fill=Status)) +
    theme_light(base_size=16) +
    scale_fill_brewer(palette="Paired") +
    coord_flip()
    #  +
    # scale_x_discrete(guide=guide_axis(n.dodge=2))

pg1=pg0 + ggtitle("Fragments in MACS Peaks") +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))

pg2=pg0 + geom_bar(stat="identity",position="fill") +
    ggtitle("Fraction of Fragments in MACS Peaks") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab("Percentage")

dd=read_tsv("out/macs/peaks_raw_fcCounts.txt",comment="#",show_col_types=FALSE)
colnames(dd)=c(colnames(dd)[1:6],pbam_to_sample(colnames(dd)[-(1:6)],samples))

#
# PCA of the samples outside EXC groups, columns renamed to SampleID
#
pca_manifest <- manifest |> filter(!str_detect(Group, "^EXC"))

pca_pages <- list()
if (nrow(pca_manifest) < 2) {
    cat("\n    No PCA:", nrow(pca_manifest),
        "sample(s) outside EXC groups; a PCA needs at least 2\n\n")
} else {

    d <- dd |>
        select(PeakNo = Geneid,
               all_of(set_names(pca_manifest$MapID, pca_manifest$SampleID))) |>
        data.frame(check.names = FALSE) |>
        column_to_rownames("PeakNo")

    group <- fct_inorder(pca_manifest$Group)
    y <- DGEList(counts=d,group=group)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)

    pr=prcomp(cpm(y,log=T),scale=F)
    dp <- pr$rotation |>
        data.frame() |>
        rownames_to_column("SampleID") |>
        left_join(pca_manifest, by = "SampleID") |>
        mutate(Group = fct_inorder(Group))

    #
    # cols25 has 25 colors; more groups (a default manifest of a large
    # run) get evenly spaced hues
    #
    nGroups <- nlevels(dp$Group)
    groupColors <- if (nGroups <= 25) cols25(nGroups) else hue_pal()(nGroups)

    pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4) + scale_color_manual(values=groupColors)
    pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_manual(values=groupColors) + geom_label_repel(color="black",max.overlaps=nrow(dp),force=10,force_pull=.1,max.iter=100000)

    pca_pages <- list(pp1, pp2)
}

dir_create(QC_DIR)
pfile <- qc_report_file("pdf", RUNTAG)
pdf(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
walk(pca_pages, print)
dev.off()

cat("Wrote", pfile, "\n")
