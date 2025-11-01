args <- commandArgs(trailing = TRUE)
if (len(args) < 3) {
  cat("
    Usage: analyzeATAC.R GENOME SampleManifest.csv Comparisons.csv [RUNTAG]

        Comparisons.csv (NO Column names)

            aNSC_loxp15,aNSC_p53

            Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15


")
  quit()
}

#
# Example inputs:
#   SampleManifest.csv
#
#        SampleID,Group,MapID
#        aNSC_loxp15_1,aNSC_loxp15,s_aNSC_loxp15_1
#        aNSC_loxp15_2,aNSC_loxp15,s_aNSC_loxp15_2
#        aNSC_p53_1,aNSC_p53,s_aNSC_p53_1
#
#   Comparisons.csv
#
#        aNSC_loxp15,aNSC_p53
#
#   Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15
#

#' Fix sample names based on naming convention
#'
#' @param ss Character vector of sample names
#' @return Character vector of fixed sample names
fix_sample_names <- function(ss) {
  if (grepl("___MD", ss[1])) {
    fix_sample_names_pemap(ss)
  } else {
    fix_sample_names_bic(ss)
  }
}

#' Fix sample names using BIC convention
#'
#' @param ss Character vector of sample names
#' @return Character vector of fixed sample names
fix_sample_names_bic <- function(ss) {
  ss %>%
    gsub("_postProcess.*", "", .) %>%
    gsub(".*_s_", "s_", .) %>%
    (\(x) sampRename[x])() %>%
    unname()
}

#' Fix sample names using PEmap convention
#'
#' @param ss Character vector of sample names
#' @return Character vector of fixed sample names
fix_sample_names_pemap <- function(ss) {
  ss %>%
    basename() %>%
    gsub("___.*", "", .) %>%
    (\(x) sampRename[x])() %>%
    unname()
}

#' Create reverse log transformation for ggplot2
#'
#' @param base Numeric base for logarithm (default: e)
#' @return A transformation object for ggplot2 scales
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(
    paste0("reverselog-", format(base)),
    trans,
    inv,
    log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

#' Create PNG file with Cairo graphics device
#'
#' @param filename Character path to output file
#' @param width Numeric width in inches (default: 14)
#' @param height Numeric height in inches (default: 8.5)
#' @param pointsize Numeric point size for text (default: 12)
#' @param res Numeric resolution in DPI (default: 150)
png_cairo <- function(filename, width = 14, height = 8.5, pointsize = 12, res = 150) {
  png(
    filename,
    type = "cairo",
    units = "in",
    width = width,
    height = height,
    pointsize = pointsize,
    res = res
  )
}


#' Merge PNG files into a single PDF
#'
#' @param file_spec Character file specification with printf-style formatting
merge_pngs <- function(file_spec) {
  file_re <- gsub("_%\\d+d", ".*", file_spec)
  pdf_file <- gsub("_%\\d+d.*", ".pdf", file_spec)
  system2(
    "convert",
    c(sort(dir_ls(regex = file_re)), pdf_file),
    stderr = cc("stderr", "mergePNGs", "convert", DATE())
  )
}

# Load required libraries
library(ChIPseeker)
library(AnnotationDbi)
library(patchwork)
library(scales)
library(edgeR)
library(ggrepel)
library(ggsci)
library(tidyverse)
library(fs)
library(openxlsx)
library(ggrastr)

# Read feature counts summary
ds <- read_tsv("peaks_raw_fcCounts.txt.summary")

# Parse command line arguments
GENOME <- args[1]
MANIFEST_FILE <- args[2]
COMPARISON_FILE <- args[3]
RUNTAG <- if (len(args) == 4) args[4] else ""

# Set genome-specific annotation databases
if (GENOME == "human") {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  annoDb <- "org.Hs.eg.db"
} else if (GENOME == "mouse") {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  annoDb <- "org.Mm.eg.db"
} else {
  cat("\n\tUnknown GENOME", GENOME, "\n")
  cat("\tValid genomes: human, mouse\n\n")
  quit()
}

# Load sample manifest and create renaming lookup
manifest <- read_csv(MANIFEST_FILE) |> arrange(SampleID)

sampRename <- manifest$SampleID
names(sampRename) <- manifest$MapID

# Process feature counts summary
ds <- ds |>
  gather(Sample, Count, -Status) |>
  mutate(Sample = fix_sample_names(Sample)) |>
  mutate(Status = gsub("_.*", "", Status)) |>
  group_by(Sample, Status) |>
  summarize(Counts = sum(Count), .groups = "drop") |>
  mutate(Status = ifelse(Status == "Assigned", "InPeaks", "Outside"))

# Load raw feature counts
dd <- read_tsv("peaks_raw_fcCounts.txt", comment = "#")

# Extract peak annotation
peak.annote <- dd |> select(PeakNo = Geneid, Chr, Start, End, Strand, Length)

# Remove excluded samples
manifest <- manifest |> filter(!grepl("^EXC", Group))

# Create count matrix
d <- dd |>
  select(PeakNo = Geneid, matches(".bam$")) |>
  data.frame(check.names = FALSE) |>
  column_to_rownames("PeakNo")

colnames(d) <- fix_sample_names(colnames(d))

# Validate count matrix columns match manifest
if (!all(colnames(d) %in% manifest$MapID)) {
  cat("\nERROR in creation of count matrix 'd'\n")
  cat("LINE-141\n\n")
  rlang::abort("ERROR")
}

# Prepare count matrix for edgeR analysis
d <- d[, manifest$SampleID]
group <- factor(manifest$Group)
y <- DGEList(counts = d, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Perform PCA for quality control
pr <- prcomp(cpm(y, log = TRUE), scale = FALSE)
dp <- pr$rotation |>
  data.frame() |>
  rownames_to_column("SampleID") |>
  left_join(manifest, by = "SampleID")

# Define color palette
colors1 <- c(pal_uchicago("default")(9), pals::cols25())

# Create PCA plots
pp1 <- ggplot(dp, aes(PC1, PC2, color = Group, label = SampleID)) +
  theme_light(base_size = 16) +
  geom_point(size = 4, alpha = 0.6) +
  scale_color_manual(values = colors1)

pp2 <- ggplot(dp, aes(PC1, PC2, color = Group, label = SampleID)) +
  theme_light(base_size = 16) +
  geom_point(size = 2) +
  scale_color_manual(values = colors1) +
  geom_label_repel(
    color = "black",
    max.overlaps = Inf,
    min.segment.length = 0,
    size = 3,
    force = 16,
    max.time = 10,
    max.iter = 100000
  )

# Set up design matrix and estimate dispersion
design <- model.matrix(~0 + group)
y <- estimateDisp(y, design)

#' Perform differential peak analysis using edgeR LRT
#'
#' @param y DGEList object containing counts
#' @param design Design matrix
#' @param contrast Contrast vector
#' @param fdr_cut FDR cutoff for significance (default: 0.05)
#' @return List containing results table, comparison name, and plots
do_qlf_stats <- function(y, design, contrast, fdr_cut = 0.05) {

  cat("Num peaks =", nrow(y), "\n")
  cat("\nRefilter peaks to just those in comp\n\n")

  # Extract samples involved in this comparison
  comp_samps <- which(rowSums(design[, contrast != 0]) > 0) |> unname()

  # Subset data to comparison samples and re-filter
  g_comp <- factor(group[comp_samps])
  y_comp <- getCounts(y)[, comp_samps]
  design_comp <- model.matrix(~0 + g_comp)

  y_comp <- DGEList(counts = y_comp, group = g_comp)
  keep <- filterByExpr(y_comp, design = design_comp)
  y_comp <- y_comp[keep, , keep.lib.sizes = FALSE]
  y_comp <- calcNormFactors(y_comp)
  y_comp <- estimateDisp(y_comp, design_comp)

  # Perform likelihood ratio test
  contrast_comp <- contrast[contrast != 0]
  fit <- glmFit(y_comp, design_comp)
  lrt <- glmLRT(fit, contrast = contrast_comp)
  model <- lrt

  # Create comparison tag
  cp <- model$comparison
  comp_tag <- paste0(sort(strsplit(cp, " ")[[1]], decreasing = TRUE), collapse = "") %>%
    gsub("-1\\*", "-", .) %>%
    gsub("^1\\*", "", .) %>%
    gsub("gComp", "", .)

  # Process results table
  tt <- model$table |>
    data.frame() |>
    rownames_to_column("PeakNo") |>
    tibble() |>
    mutate(FDR = p.adjust(PValue)) |>
    arrange(FDR, PValue) |>
    mutate(PValue.mod = ifelse(PValue < .Machine$double.eps^2, .Machine$double.eps^2, PValue))

  # Apply IHW (Independent Hypothesis Weighting) for multiple testing correction
  library(IHW)
  ihw_res <- ihw(PValue ~ Length, dat = tt |> left_join(peak.annote, by = "PeakNo"), alpha = 0.25)
  cat("ihw rejections =", rejections(ihw_res), "\n")

  # Create MA plot
  max_logfc <- max(abs(tt$logFC))

  p_ma <- ggplot(arrange(tt, desc(PValue)), aes(logCPM, logFC, color = FDR < fdr_cut)) +
    theme_light(base_size = 16) +
    rasterize(geom_point()) +
    scale_color_manual(values = c("#7f7f7f33", "#e31a1c")) +
    ggtitle(comp_tag) +
    scale_y_continuous(limits = c(-1, 1) * max_logfc)

  # Create volcano plot (PValue.mod is clipped to make plot readable)
  p_vc <- ggplot(arrange(tt, desc(PValue)), aes(logFC, PValue.mod, color = FDR < fdr_cut)) +
    theme_light(base_size = 16) +
    rasterize(geom_point()) +
    scale_y_continuous(trans = reverselog_trans(10)) +
    scale_x_continuous(limits = c(-1, 1) * max_logfc) +
    scale_color_manual(values = c("#7f7f7f33", "#e31a1c")) +
    ggtitle(comp_tag)

  # Filter for significant peaks
  ans <- tt |>
    filter(FDR < fdr_cut) |>
    left_join(peak.annote, by = "PeakNo") |>
    select(-matches("^F$|^LR")) |>
    filter(Chr != "MT")

  # Prepare significant peaks for annotation
  sig_peaks <- ans |>
    select(V1 = Chr, V2 = Start, V3 = End, V4 = PeakNo, V5 = PValue) |>
    mutate(V5 = -10 * log10(V5))

  # Annotate significant peaks with gene information
  if (nrow(ans) > 0) {
    # Add "chr" prefix if missing
    if (substr(ans$Chr[1], 1, 3) != "chr") {
      sig_peaks <- sig_peaks |> mutate(V1 = paste0("chr", V1)) |> data.frame()
    } else {
      sig_peaks <- sig_peaks |> data.frame()
    }

    # Convert to GRanges and annotate
    sig_peaks <- ChIPseeker:::peakDF2GRanges(sig_peaks)
    aa <- annotatePeak(sig_peaks, TxDb = txdb, tssRegion = c(-5000, 5000), annoDb = annoDb)

    # Extract gene annotations
    peak_annote_genes <- as.data.frame(aa) |>
      tibble() |>
      dplyr::select(PeakNo = V4, SYMBOL, GENENAME, annotation, distanceToTSS)

    tbl <- ans |> left_join(peak_annote_genes, by = "PeakNo")
  } else {
    cat("\n   No significant peaks", comp_tag, "at FDR", fdr_cut, "\n\n")
    tbl <- ans
  }

  list(tbl = tbl, comparison = comp_tag, p.ma = p_ma, p.vc = p_vc, model = model)
}

# Extract project number from working directory
cwd <- strsplit(getwd(), "/")[[1]]
proj_no <- grep("^Proj_|^B-\\d+", cwd, value = TRUE)

# Process comparisons
g_names <- gsub("group", "", colnames(design))
comps <- read_csv(COMPARISON_FILE, col_names = FALSE)

res <- list()

for (ci in transpose(comps)) {
  contrast <- as.numeric(g_names == ci$X2) - as.numeric(g_names == ci$X1)
  cat("========================================================\n")
  cat(str(ci))

  # Validate contrast
  if (!(sum(contrast != 0) > 0 & sum(contrast) == 0)) {
    cat("\n\n")
    cat("   Invalid contrasts")
    cat("\n\n")
    stop("FATAL")
  }

  res[[len(res) + 1]] <- do_qlf_stats(y, design, contrast)
}
cat("========================================================\n")
cat("========================================================\n\n")

# Compile results and summary statistics
tbls <- map(res, "tbl")
names(tbls) <- map(res, "comparison") |> unlist() |> substr(1, 31)
stats <- map(tbls, nrow) |> bind_rows() |> gather(Comparison, NumSig)

# Write results to Excel file
write.xlsx(
  c(list(Summary = stats), tbls),
  cc(proj_no, RUNTAG, "DiffPeaksEdgeR_V3.xlsx")
)

# Generate PDF output with all plots
pfile <- cc(proj_no, RUNTAG, "DiffPeaks_V3.pdf")
pdf(pfile, width = 11, height = 8.5)

print(pp1)
print(pp2)
for (ii in seq(len(res))) {
  print(res[[ii]]$p.ma)
  print(res[[ii]]$p.vc)
}
dev.off()


