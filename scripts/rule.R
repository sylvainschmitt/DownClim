# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
folderout <- snakemake@output[[1]]

# libraries
suppressMessages(library(tidyverse)) 

# code
dir.create(folderout) # dummy
