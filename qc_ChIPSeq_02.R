#
# qc_ChIPSeq_02.R MANIFEST [RUNTAG]
#
# QC stage 2, run by postPipe.sh from the analysis directory when
# manifest.txt exists. MANIFEST has no header and two columns,
# SampleID<TAB>Group; every SampleID must be a sample of the run (its
# @RG SM tag) exactly. Groups starting with EXC are left out of the PCA.
#
# Writes qcChIPSeq_<Proj>_<RUNTAG>_.pdf: fragments in and outside the
# merged MACS peaks per sample, and a PCA of the peak counts.
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

if(length(args)<1) {
    fatal("usage: qc_ChIPSeq_02.R MANIFEST [RUNTAG]")
}
if(length(args)==2) {
    RUNTAG=args[2]
} else {
    RUNTAG=""
}

samples=read_samples()

manifest=read_tsv(args[1],col_names=c("SampleID","Group"),
                  col_types=cols(.default="c"),show_col_types=FALSE)

notSamples=setdiff(manifest$SampleID,samples$SAMPLE)
if(length(notSamples)>0) {
    fatal("manifest SampleIDs that are not samples of this run:",
          str_c(notSamples,collapse=" "),
          "; samples:",str_c(samples$SAMPLE,collapse=" "))
}

ds=read_tsv("out/macs/peaks_raw_fcCounts.txt.summary",show_col_types=FALSE)
#
# Columns are the postProcess BAM paths; samples.tsv maps them to sample
# names exactly.
#
colnames(ds)=c("Status",pbam_to_sample(colnames(ds)[-1],samples))

ds=ds %>%
    gather(Sample,Count,-Status) %>%
    mutate(Status=str_remove(Status,"_.*")) %>%
    group_by(Sample,Status) %>%
    summarize(Counts=sum(Count),.groups="drop") %>%
    mutate(Status=ifelse(Status=="Assigned","InPeaks","Outside"))

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

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

d=dd %>%
    select(PeakNo=Geneid,all_of(manifest$SampleID)) %>%
    data.frame(check.names=F) %>%
    column_to_rownames("PeakNo")

#
# Remove excluded points
#

manifest=manifest %>% filter(!str_detect(Group,"^EXC"))

d=d[,manifest$SampleID]
group=factor(manifest$Group)
y <- DGEList(counts=d,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

pr=prcomp(cpm(y,log=T),scale=F)
dp=pr$rotation %>% data.frame %>% rownames_to_column("SampleID") %>% left_join(manifest,by="SampleID")

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4) + scale_color_manual(values=cols25())
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_manual(values=cols25()) + geom_label_repel(color="black",max.overlaps=nrow(manifest),force=10,force_pull=.1,max.iter=100000)

projNo=get_project_number()

pfile=str_c("qcChIPSeq",projNo,RUNTAG,".pdf",sep="_")
pdf(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
print(pp1)
print(pp2)
dev.off()
