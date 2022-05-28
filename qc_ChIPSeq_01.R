require(tidyverse)
require(fs)
require(openxlsx)

peakFiles=dir_ls("out/macs",recur=T,regex="peaks.xls$")

qCut=0.05

#
# Get stats and write sig peak file
#

getQCAndWriteSigPeakFile<-function(ff,qCut=0.05) {

    pbase=basename(ff)

    dp=read_tsv(ff,comment="#",col_types=cols(chr=col_character()))
    dp.sig=dp %>% arrange(desc(`-log10(qvalue)`)) %>% filter(`-log10(qvalue)`>-log10(qCut))

    write.xlsx(
        dp.sig,
        gsub("_peaks.xls",paste0("___sigPeaks__FDR_",qCut,".xlsx"),pbase)
        )

    png(gsub("_peaks.xls",paste0("___Volcano__FDR_",qCut,".png"),pbase),type="cairo",
            units="in", width=8, height=8, pointsize=12, res=150)

    sampName=basename(dirname(ff))
    numPeaks=nrow(dp)
    numSigPeaks=nrow(dp.sig)

    kk=dp$`-log10(qvalue)`>-log10(qCut)

    plot(dp$fold_enrichment,10^(-dp$`-log10(pvalue)`),log="xy",
        xlim=c(1,max(dp$fold_enrichment)),
        xlab="FoldEnrichment",ylab="PValue",
        main=sampName)
    points(dp$fold_enrichment[kk],10^(-dp$`-log10(pvalue)`)[kk],pch=19,col="darkred")
    dev.off()

    tibble(
        Sample=sampName,
        TotalPeaks=numPeaks,
        FDR=qCut,
        SigPeaks=numSigPeaks,
        PCT=SigPeaks/TotalPeaks
        )

}

xx=map(peakFiles,getQCAndWriteSigPeakFile)
stats=bind_rows(xx)

projNo=unique(gsub("_s_.*","",basename(dirname(peakFiles))))

write.xlsx(stats,cc("qcChIPSeq",projNo,".xlsx"))
