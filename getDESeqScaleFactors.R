#
# getDESeqScaleFactors.R COUNTS
#
# DESEQ stage: write scaleFactorsDESeq2.csv beside COUNTS (the
# featureCounts matrix), holding the inverse of DESeq2's sizeFactors.
# The inversion is deliberate (NOTES.md).
#

cat("## START getDESeqScaleFactors.R", format(Sys.time()), "\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    cat("\n    FATAL ERROR: usage: getDESeqScaleFactors.R COUNTS\n\n", file = stderr())
    quit(save = "no", status = 1)
}

suppressPackageStartupMessages(library(DESeq2))

dd=read.delim(args[1],comment="#")
ds=as.matrix(dd[,7:ncol(dd)])
rownames(ds)=dd$Geneid
cond=factor(colnames(ds))
dds=DESeqDataSetFromMatrix(ds,colData=data.frame(Sample=cond),design=~Sample)
dds=estimateSizeFactors(dds)
sfm=data.frame(scaleFactor=1/sizeFactors(dds))
write.csv(sfm,file.path(dirname(args[1]),"scaleFactorsDESeq2.csv"))

cat("## END getDESeqScaleFactors.R", format(Sys.time()), "\n")
