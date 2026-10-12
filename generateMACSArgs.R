#
# generateMACSArgs.R PAIRS SAMPLES_TSV
#
# Resolve the pairing file (control in column 1, target in column 2, one
# pair per line; 'na' in column 1, in any case, for no control) against
# the run's samples (out/pipeline_info/samples.tsv). Prints one line per
# pair for out/macsPairs.txt:
#
#   out/<target>_postProcess.clean.bed.gz out/<control>_postProcess.clean.bed.gz
#   out/<target>_postProcess.clean.bed.gz _NoCTRL
#
# Names must equal a sample name (the @RG SM tag) exactly. A name that
# matches no sample, a missing BED, a malformed or repeated line is fatal,
# with the line of the pairing file.
#

suppressPackageStartupMessages({
    library(stringr)
})
SDIR <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |> str_remove("^--file=") |> dirname()
source(file.path(SDIR, "tools.R"))

cat("## START generateMACSArgs.R", format(Sys.time()), "\n", file = stderr())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    fatal("usage: generateMACSArgs.R PAIRS SAMPLES_TSV")
}
pairingFile <- args[1]
samplesFile <- args[2]

if (!file_exists(pairingFile)) {
    fatal(str_glue("no pairing file [{pairingFile}]"))
}

samples <- read_samples(samplesFile)

lines <- read_lines(pairingFile)
lineNo <- seq_along(lines)
keep <- str_detect(lines, "\\S")
lines <- lines[keep]
lineNo <- lineNo[keep]

if (length(lines) == 0) {
    fatal(str_glue("pairing file [{pairingFile}] is empty"))
}

bed_of <- function(sample) {
    str_glue("out/{sample}_postProcess.clean.bed.gz")
}

out <- character(0)

for (ii in seq_along(lines)) {

    where <- str_glue("{pairingFile} line {lineNo[ii]} [{lines[ii]}]")
    fields <- str_split_1(str_trim(lines[ii]), "\\s+")

    if (length(fields) != 2) {
        fatal(str_glue("{where}: need two columns, CONTROL TARGET"))
    }

    control <- fields[1]
    target <- fields[2]

    if (!(target %in% samples$SAMPLE)) {
        fatal(str_glue("{where}: target [{target}] is not a sample of this run; samples: {str_c(samples$SAMPLE, collapse = ' ')}"))
    }
    tbed <- bed_of(target)
    if (!file_exists(tbed)) {
        fatal(str_glue("{where}: no BED for target [{tbed}]"))
    }

    if (str_to_lower(control) == "na") {
        cbed <- "_NoCTRL"
    } else {
        if (!(control %in% samples$SAMPLE)) {
            fatal(str_glue("{where}: control [{control}] is not a sample of this run; use 'na' for no control; samples: {str_c(samples$SAMPLE, collapse = ' ')}"))
        }
        if (control == target) {
            fatal(str_glue("{where}: target is its own control"))
        }
        cbed <- bed_of(control)
        if (!file_exists(cbed)) {
            fatal(str_glue("{where}: no BED for control [{cbed}]"))
        }
    }

    pair <- str_c(tbed, cbed, sep = " ")
    if (pair %in% out) {
        fatal(str_glue("{where}: pair given twice"))
    }
    out <- c(out, pair)
}

cat(out, sep = "\n")

cat("## END generateMACSArgs.R", format(Sys.time()), "\n", file = stderr())
