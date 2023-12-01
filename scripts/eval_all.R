# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
filesin <- snakemake@input[[1]]
fileout <- snakemake@output[[1]]

# libraries
library(tidyverse)

# code
lapply(filesin, read_tsv) %>% 
  bind_rows() %>% 
  write_tsv(fileout)
