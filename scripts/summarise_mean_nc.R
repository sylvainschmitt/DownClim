# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
filein <- snakemake@input[[1]]
fileout <- snakemake@output[[1]]
years <- snakemake@params$years

# test
# filein <- "results/cordex/cropped/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc"
# years <- 2000:2010

# libraries
library(tidyverse)
library(terra)

# code
r <- rast(filein)
subset(r, year(time(r)) %in% years) %>% 
  tapp(months, mean) %>% 
  writeCDF(fileout)
