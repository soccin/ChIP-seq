suppressPackageStartupMessages({
    require(tidyverse)
})

source("get_summit_windows.R")

args=commandArgs(trailing=T)
summitFile=args[1]

bedCols=cols(
  X1 = col_character(),
  X2 = col_double(),
  X3 = col_double(),
  X4 = col_character(),
  X5 = col_double()
)

MAX_PEAKS=5000

bed=read_tsv(summitFile,col_names=F,col_types=bedCols)

xx=bed |> arrange(desc(X5)) %>% slice(1:min(5000,nrow(.))) |> transpose() |> map(get_summit_windows,W=200,numBins=50) |> bind_rows()
write_tsv(xx,gsub(".bed$","__Profile.bed",basename(summitFile)),col_names=F)



