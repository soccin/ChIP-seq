#
# getMedianFragmentLengthFromPredictDFile.R LOG [LOG ...]
#
# Print the median (floored) of the fragment lengths that macs2 predictd
# wrote into the LOGs (out/profiles/*.log). Samples without an estimate
# are listed on stderr. Fatal if no LOG has one.
#

suppressPackageStartupMessages({
    library(stringr)
    library(purrr)
})
SDIR <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |> str_remove("^--file=") |> dirname()
source(file.path(SDIR, "tools.R"))

cat("## START getMedianFragmentLengthFromPredictDFile.R", format(Sys.time()), "\n", file = stderr())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    fatal("usage: getMedianFragmentLengthFromPredictDFile.R LOG [LOG ...]")
}

fragLen <- map_dbl(args, \(f) {
    read_lines(f) |>
        str_extract("# predicted fragment length is (\\d+) bps", group = 1) |>
        discard(is.na) |>
        as.numeric() |>
        first(default = NA_real_)
})

walk2(args, fragLen, \(f, n) cat("   ", f, "fragment length", n, "\n", file = stderr()))

if (any(is.na(fragLen))) {
    cat("    no estimate in:", args[is.na(fragLen)], "\n", file = stderr())
}

if (all(is.na(fragLen))) {
    fatal("no fragment length estimate in any predictd log")
}

med.FragLen <- floor(median(fragLen, na.rm = TRUE))
cat(med.FragLen, "\n")

cat("## END getMedianFragmentLengthFromPredictDFile.R", format(Sys.time()), "\n", file = stderr())
