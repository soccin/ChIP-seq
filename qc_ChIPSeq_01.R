suppressPackageStartupMessages({
require(tidyverse)
require(fs)
require(openxlsx)
})

source("ChIP-seq/tools.R")

peakFiles=dir_ls("out/macs",recur=T,regex="peaks.xls$")

qCut=0.05

#
# Get stats and write sig peak file
#

getQCAndWriteSigPeakFile<-function(ff,qCut=0.05) {

    pbase=basename(ff)

    cat("File =",ff,"\n")

    dp=read_tsv(ff,comment="#",col_types=cols(chr=col_character()))
    dp.sig=dp %>% arrange(desc(`-log10(qvalue)`)) %>% filter(`-log10(qvalue)`>-log10(qCut))

    sampName=basename(dirname(ff))
    numPeaks=nrow(dp)
    numSigPeaks=nrow(dp.sig)

    if(numPeaks>0) {

        if(numSigPeaks>0) {
            write.xlsx(
                dp.sig,
                gsub("_peaks.xls",paste0("___sigPeaks__FDR_",qCut,".xlsx"),pbase)
                )
        }

        png(gsub("_peaks.xls",paste0("___Volcano__FDR_",qCut,".png"),pbase),type="cairo",
                units="in", width=8, height=8, pointsize=12, res=150)

        numPeaks=nrow(dp)
        numSigPeaks=nrow(dp.sig)

        kk=dp$`-log10(qvalue)`>-log10(qCut)

        Qscore=pmax(10^(-dp$`-log10(pvalue)`),.Machine$double.xmin)

        plot(dp$fold_enrichment,Qscore,log="xy",
            xlim=c(1,max(dp$fold_enrichment)),
            xlab="FoldEnrichment",ylab="PValue",
            main=sampName)
        points(dp$fold_enrichment[kk],Qscore[kk],pch=19,col="darkred")
        dev.off()

    }

    tibble(
        Sample=sampName,
        TotalPeaks=numPeaks,
        FDR=qCut,
        SigPeaks=numSigPeaks,
        PCT=ifelse(numPeaks>0,SigPeaks/TotalPeaks,0)
        )


}

xx=map(peakFiles,getQCAndWriteSigPeakFile)
stats=bind_rows(xx)

extract_common_prefix<-function(x) {

    # sort the vector
    x<-sort(x)
    # split the first and last element by character
    d_x<-strsplit(x[c(1,length(x))],"")
    minLen=min(map_vec(d_x,len))
    d_x=map(d_x,~.[1:minLen])
    # search for the first not common element and so, get the last matching one
    der_com<-match(FALSE,do.call("==",d_x))-1
    # if there is no matching element, return an empty vector, else return the common part
    ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))

}

projNo=get_project_number()

write.xlsx(stats,cc("qcChIPSeq",projNo,".xlsx"))

