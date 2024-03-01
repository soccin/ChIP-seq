suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(purrr))
suppressPackageStartupMessages(require(readr))
args=commandArgs(trailing=T)
pairingFile=args[1]
inputs=args[-1]

pairs=read_tsv(pairingFile,col_names=F,col_type=cols())

for(ii in seq(nrow(pairs))) {

    target=grep(paste0(pairs$X2[ii],"_.*postProcess"),inputs,value=T)

    if(len(target)!=1) {
        cat("\n\nCan not find target sample",pairs$X2[ii],"\n",file=stderr())
        cat("\n\n",file=stderr())
        rlang::abort("FATAL::ERROR")
    }

    if(!(pairs$X1[ii] %in% c("NA","na","_na","_NA"))) {

        control=grep(paste0(pairs$X1[ii],"_postProcess"),inputs,value=T)

        if(len(control)==0) {
            control="_NoCTRL"
        }

    } else {
        control="_NoCTRL"
    }

    if(len(target)>1 | len(control)>1) {
        cat("\n\tFATAL ERROR: degenerate target/control\n\n")
        cat(target)
        cat("\n")
        cat(control)
        cat("\n\n")
        stop("FATAL ERROR")
    }

    cat(target,control,"\n")

}
