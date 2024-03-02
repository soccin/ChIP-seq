args=commandArgs(trailing=T)
if(len(args)<1) {
     cat("\n   Usage: analyzeATAC.R SampleManifest.txt [RUNTAG]\n\n")
     quit()
}
if(len(args)==2) {
    RUNTAG=args[2]
} else {
    RUNTAG=""
}

suppressPackageStartupMessages({
require(tidyverse)
require(scales)
require(patchwork)
require(edgeR)
require(ggrepel)
require(pals)
})

source("ChIP-seq/tools.R")

manifest=read_tsv(args[1],col_names=c("SampleID","Group"))

# require(fs)
# require(openxlsx)

ds=read_tsv("out/macs/peaks_raw_fcCounts.txt.summary")
#
# Clean up the samples names. Work with BIC pipeline
# and PEMapper filename conventions:
# - BIC: Proj_15426_s_IgGnegcontrol1_postProcess.bam
# - PEMapper: s_Scramble3-IgG___MD_postProcess.bam
#

fix_colnames<-function(dx) {
    basename(colnames(dx)) %>% 
        gsub(".*_s_","s_",.) %>% 
        gsub("_postProcess.bam$","",.) %>% 
        gsub("___MD$","",.)
}
colnames(ds)=fix_colnames(ds)

ds=ds %>%
    gather(Sample,Count,-Status) %>%
    mutate(Status=gsub("_.*","",Status)) %>%
    group_by(Sample,Status) %>%
    summarize(Counts=sum(Count)) %>%
    mutate(Status=ifelse(Status=="Assigned","InPeaks","Outside"))

pg0=ggplot(ds,aes(Sample,Counts,fill=Status)) +
    theme_light(base_size=16) +
    scale_fill_brewer(palette="Paired") +
    coord_flip()
    #  +
    # scale_x_discrete(guide=guide_axis(n.dodge=2))

pg1=pg0 + ggtitle("Mapped Reads in MACS Peaks") +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))

pg2=pg0 + geom_bar(stat="identity",position="fill") +
    ggtitle("Fraction Reads in MACS Peaks") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab("Percentage")

dd=read_tsv("out/macs/peaks_raw_fcCounts.txt",comment="#")
colnames(dd) = fix_colnames(dd)

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

d=dd %>%
    select(PeakNo=Geneid,all_of(manifest$SampleID)) %>%
    data.frame(check.names=F) %>%
    column_to_rownames("PeakNo")

#
# Remove excluded points
#

manifest=manifest %>% filter(!grepl("^EXC",Group))

d=d[,manifest$SampleID]
group=factor(manifest$Group)
y <- DGEList(counts=d,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

pr=prcomp(cpm(y,log=T),scale=F)
dp=pr$rotation %>% data.frame %>% rownames_to_column("SampleID") %>% left_join(manifest)

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4) + scale_color_manual(values=cols25())
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_manual(values=cols25()) + geom_label_repel(color="black",max.overlaps=nrow(manifest),force=10,force_pull=.1,max.iter=100000)

projNo=get_project_number()

pfile=cc("qcChIPSeq",projNo,RUNTAG,".pdf")
pdf(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
print(pp1)
print(pp2)
dev.off()


