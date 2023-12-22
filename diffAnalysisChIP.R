require(tidyverse)
require(edgeR)

projNo=grep("Proj_",strsplit(fs::path_real("pipeline"),"/")[[1]],value=T)

dd=read_tsv("out/macs/peaks_raw_fcCounts.txt",comment="#") %>%
    rename_all(~basename(.)%>%gsub(".*_s_","s_",.)%>%gsub("_postProcess.*","",.))

key=read_csv("groups.csv")
comps=read_tsv("comparisons.txt",col_names=F) %>% transpose

qCut=0.05

annote=read_tsv("annote/macsPeaksMerged.bed_HOMERAnnote.txt") %>%
    select(-Chr,-Start,-End,-Strand,-`Peak Score`,-`Focus Ratio/Region Size`)
colnames(annote)[1]="PeakID"

for(ci in comps) {

    cat("Comparison",paste0(ci$X1,"-",ci$X2),"\n")

    grps=unname(unlist(ci))
    key0=key %>% filter(Group %in% grps)

    ds=dd %>% select(Geneid,all_of(key0$Sample)) %>% column_to_rownames("Geneid")

    groups=key0$Group
    names(groups)=key0$Sample
    groups=groups[colnames(ds)]
    groups=factor(groups,levels=c(ci$X1,ci$X2))

    design=model.matrix(~0+groups)

    y <- DGEList(counts=ds,group=groups)

    keep <- filterByExpr(y,design,min.count=8)
    y <- y[keep, , keep.lib.sizes=FALSE]

    y <- calcNormFactors(y)

    y <- estimateDisp(y,design,robust=T)

    fit <- glmQLFit(y, design)

    qlf <- glmQLFTest(fit, contrast=c(1,-1))

    nSig=sum(topTags(qlf,n=nrow(y))$table$FDR<qCut)

    tbl=topTags(qlf,n=nSig)$table %>%
        rownames_to_column("PeakID") %>%
        tibble %>%
        select(PeakID,PValue,FDR,logFC) %>%
        mutate(FC=ifelse(logFC>0,2^logFC,-(2^(-logFC))))

    medianCounts=cpm(y) %>%
        as.data.frame %>%
        rownames_to_column("PeakID") %>%
        tibble %>%
        gather(Sample,Count,-PeakID) %>%
        left_join(key0) %>%
        mutate(Group=factor(Group,levels=c(ci$X1,ci$X2))) %>%
        group_by(PeakID,Group) %>%
        summarize(medianCount=median(Count)) %>%
        spread(Group,medianCount)

    tbl=left_join(tbl,medianCounts) %>%
        left_join(dd %>% select(PeakID=Geneid,Chr,Start,End)) %>%
        left_join(annote)

    openxlsx::write.xlsx(tbl,cc(projNo,"_diffAnalysis_",paste0(ci$X1,"-",ci$X2),".xlsx"))

}

