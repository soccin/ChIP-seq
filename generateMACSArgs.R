suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(purrr))
suppressPackageStartupMessages(require(readr))
args=commandArgs(trailing=T)
pairingFile=args[1]
inputs=args[-1]

pairs=read_tsv(pairingFile,col_names=F,col_type=cols())

for(ii in seq(nrow(pairs))) {

    target=grep(pairs$X2[ii],inputs,value=T)
    control=grep(pairs$X1[ii],inputs,value=T)

    if(len(control)==0) {
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
