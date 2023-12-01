# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
proj_file <- snakemake@input[[1]]
base_file <- snakemake@input[[2]]
fileout <- snakemake@output[[1]]
var <- snakemake@params$var

# test
# proj_file <- "results/cordex/proj/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc"
# base_file <- "results/cordex/hist/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc"

# libraries
library(tidyverse)
library(terra)

# code
proj <- rast(proj_file)
base <- rast(base_file)
relative <- FALSE
if(var == "pr")
  relative = TRUE
if(relative)
  warning("Using relative bias.")
if(relative)
  anomaly <- (proj - base)/base
if(!relative)
  anomaly <- proj - base
writeCDF(x = anomaly, filename = fileout)
