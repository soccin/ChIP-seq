extractSampleNames<-function(ss,pno,suffix) {
    gsub(paste0(".*",pno,"_"),"",ss) %>% gsub(suffix,"",.)
}

#halt("INCLUDE")

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
})

manifest=read_tsv(args[1],col_names=c("SampleID","Group"))

# require(fs)
# require(openxlsx)

ds=read_tsv("out/macs/peaks_raw_fcCounts.txt.summary")

#
# Extract Project No
#

ci=grep(manifest$SampleID[1],colnames(ds),value=T)
projNo=basename(ci) %>% gsub(paste0("_",manifest$SampleID[1],".*"),"",.)

ds=ds %>%
    gather(Sample,Count,-Status) %>%
    mutate(Sample=extractSampleNames(Sample,projNo,"_postProcess.bam")) %>%
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

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

d=dd %>%
    select(PeakNo=Geneid,matches("Proj_.*post")) %>%
    data.frame(check.names=F) %>%
    column_to_rownames("PeakNo")

colnames(d)=extractSampleNames(colnames(d),projNo,"_postProcess.bam")

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

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4) + scale_color_brewer(palette="Paired")
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_brewer(palette="Paired") + geom_label_repel(color="black",max.overlaps=nrow(manifest),force=10,force_pull=.1,max.iter=100000)

pfile=cc("qcChIPSeq",projNo,RUNTAG,".pdf")
pdf(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
print(pp1)
print(pp2)
dev.off()


