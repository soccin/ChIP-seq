suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(purrr))
args=commandArgs(trailing=T)
med.FragLen=map(args,readLines) %>%
    map(function(x){grep("# predicted fragment length",x,value=T)}) %>%
    map(function(x){strsplit(x," ")[[1]][14]}) %>%
    unlist %>%
    as.numeric %>%
    median
cat(med.FragLen,"\n")