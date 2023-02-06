require(tidyverse)
require(ComplexHeatmap)
require(circlize)

args=commandArgs(trailing=T)
profileFile=args[1]

htitle=gsub("_summits.*","",strsplit(profileFile,"_s_")[[1]][2]) %>% gsub("____","/",.)


xx=read_tsv(profileFile,col_names=F)
mm=xx %>%
    select(-X1,-X2,-X3) %>%
    spread(X6,X7) %>%
    arrange(desc(X5)) %>%
    select(-X5) %>%
    data.frame %>%
    column_to_rownames("X4") %>%
    as.matrix

colors=colorRamp2(seq(quantile(mm,0.01),quantile(mm,.99),len=9),RColorBrewer::brewer.pal(9, "Reds"))

pdf(file=cc("pltSummitProfile",gsub("/","___",htitle),".pdf"),width=3,height=12)
Heatmap(mm,
    cluster_rows=F,cluster_columns=F,
    show_row_names=F,show_column_names=F,
    raster_by_magick=T,col=colors,
    heatmap_legend_param=list(title="",legend_height=unit(1,"npc")),
    column_title=htitle)
dev.off()

