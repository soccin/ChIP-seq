#!/usr/bin/env Rscript

# Script to extract versions of all R packages used in ChIP-seq pipeline
# Usage: Rscript getPackageVersions.R > R_requirements.txt

packages <- c(
  # CRAN packages
  'tidyverse',
  'fs',
  'openxlsx',
  'scales',
  'patchwork',
  'ggrepel',
  'ggrastr',
  'ggsci',
  'pals',
  'circlize',
  'dplyr',
  'purrr',
  'readr',
  'rlang',
  'RColorBrewer',
  # Bioconductor packages
  'edgeR',
  'DESeq2',
  'ChIPseeker',
  'AnnotationDbi',
  'ComplexHeatmap',
  'IHW',
  'org.Hs.eg.db',
  'org.Mm.eg.db',
  'TxDb.Hsapiens.UCSC.hg19.knownGene',
  'TxDb.Mmusculus.UCSC.mm10.knownGene'
)

cat("# R Package Requirements for ChIP-seq Pipeline v0.8.0\n")
cat("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("# R Version:", R.version.string, "\n\n")

cat("## Installed Packages\n\n")

installed <- c()
missing <- c()

for (pkg in sort(packages)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- as.character(packageVersion(pkg))
    cat(sprintf("%-50s %s\n", pkg, ver))
    installed <- c(installed, pkg)
  } else {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\n## Missing Packages\n\n")
  for (pkg in sort(missing)) {
    cat(sprintf("%-50s NOT_INSTALLED\n", pkg))
  }
}

cat("\n## Summary\n")
cat("Installed:", length(installed), "/", length(packages), "\n")
if (length(missing) > 0) {
  cat("Missing:", length(missing), "\n")
}
