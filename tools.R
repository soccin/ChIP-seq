#
# Shared helpers for the ChIP-seq R scripts. Each script finds this file
# from its own location and sources it:
#
#   SDIR <- commandArgs(trailingOnly = FALSE) |>
#       str_subset("^--file=") |> str_remove("^--file=") |> dirname()
#   source(file.path(SDIR, "tools.R"))     # (R/ scripts: "..")
#
# Nothing here, or in the scripts, depends on ~/.Rprofile; the pipeline
# runs them with Rscript --no-init-file.
#
# Every function that reads a run record expects the analysis directory
# as the working directory.
#

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(stringr)
    library(fs)
})

CHIP_TOOLS_DIR <- dirname(path_abs(sys.frame(1)$ofile))

RUNSTATUS_FILE <- "00.RUNSTATUS.txt"
SAMPLES_FILE <- "out/pipeline_info/samples.tsv"

#
# Print an error to stderr and exit 1.
#
fatal <- function(...) {
    cat("\n    FATAL ERROR:", ..., "\n\n", file = stderr())
    quit(save = "no", status = 1)
}

get_project_number <- function() {
    projNo <- str_split_1(getwd(), "/") |> str_subset("^Proj_") |> first()
    if (is.na(projNo)) {
        projNo <- Sys.getenv("PROJECT_NUMBER")
        cat("\n")
        cat("Can not find project number from path\n")
        cat("using env var {PROJECT_NUMBER}\n")
        cat("\n\tprojNo =", str_glue("[{projNo}]"), "\n\n")
    }
    projNo
}

#
# The genome table, lib/genomes/GENOMES.tsv.
#
genome_table <- function() {
    read_tsv(path(CHIP_TOOLS_DIR, "lib/genomes/GENOMES.tsv"),
             col_types = cols(.default = "c"), show_col_types = FALSE)
}

#
# One value of the genome table; fatal for an unknown tag or field.
#
genomeInfo <- function(tag, field) {
    tbl <- genome_table()
    if (!(tag %in% tbl$TAG)) {
        fatal(str_glue("unsupported genome [{tag}]; supported: {str_c(tbl$TAG, collapse = ' ')}"))
    }
    if (!(field %in% colnames(tbl))) {
        fatal(str_glue("unknown genome field [{field}]"))
    }
    tbl |> filter(TAG == tag) |> pull(all_of(field))
}

#
# The genome of the run in the working directory: GENOME= in
# 00.RUNSTATUS.txt, written by pipe.sh from the BAM headers. Fatal if
# there is no status file or the tag is not in the genome table.
#
run_genome <- function() {
    if (!file_exists(RUNSTATUS_FILE)) {
        fatal(str_glue("no {RUNSTATUS_FILE} in {getwd()}; run from the analysis directory of a pipe.sh run"))
    }
    genome <- read_lines(RUNSTATUS_FILE) |>
        str_subset("^GENOME=") |>
        str_remove("^GENOME=") |>
        last()
    if (is.na(genome) || genome == "") {
        fatal(str_glue("no GENOME in {RUNSTATUS_FILE}"))
    }
    genomeInfo(genome, "TAG")
}

#
# The samples of the run: out/pipeline_info/samples.tsv, written by
# pipe.sh. SAMPLE is the @RG SM tag, BAM the input BAM (absolute path),
# PBAM the postProcess BAM as it is named in the count matrix columns.
#
read_samples <- function(file = SAMPLES_FILE) {
    if (!file_exists(file)) {
        fatal(str_glue("no sample table [{file}]"))
    }
    read_tsv(file, col_types = cols(.default = "c"), show_col_types = FALSE)
}

#
# Map featureCounts column names (PBAM paths) to sample names, exactly.
# Fatal if a column is not a PBAM of the run.
#
pbam_to_sample <- function(cols, samples) {
    idx <- match(cols, samples$PBAM)
    if (anyNA(idx)) {
        fatal("count columns not in", SAMPLES_FILE, ":", str_c(cols[is.na(idx)], collapse = " "))
    }
    samples$SAMPLE[idx]
}
