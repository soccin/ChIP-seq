#
# postProcessHomer.R HOMER_ANNOTATION
#
# Write HOMER_ANNOTATION (annotatePeaks.pl output, .txt) as .xlsx, with
# the column names cleaned.
#

suppressPackageStartupMessages({
library(tidyverse)
library(openxlsx)
})

args=commandArgs(trailing=T)
if(length(args)!=1) {
    cat("\n    FATAL ERROR: usage: postProcessHomer.R HOMER_ANNOTATION\n\n",file=stderr())
    quit(save="no",status=1)
}
homerFile=args[1]

ann <- read_tsv(homerFile,show_col_types=FALSE) %>%
    rename_at(1,~str_remove(.," .*")) %>%
    rename_all(~str_replace_all(.,"[ /]","_")) %>%
    select(-Focus_Ratio_Region_Size)

write.xlsx(as.data.frame(ann),str_replace(homerFile,"\\.txt$",".xlsx"))
