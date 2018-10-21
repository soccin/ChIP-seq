library(tidyverse)
library(openxlsx)

args=commandArgs(trailing=T)
homerFile=args[1]

ann <- read_tsv(homerFile) %>%
    rename_at(1,~gsub(" .*","",.)) %>%
    rename_all(~gsub("[ /]","_",.)) %>%
    select(-Focus_Ratio_Region_Size)

write.xlsx(as.data.frame(ann),gsub(".txt",".xlsx",homerFile))
