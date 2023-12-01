# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
bias_file <- snakemake@input[[1]]
base_file <- snakemake@input[[2]]
fileout <- snakemake@output[[1]]
var <- snakemake@params$var
month <- snakemake@params$month

# test
# bias_file <- "results/cordex/proj_bias/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc"
# base_file <- "results/chelsa/hist/pr_New-Caledonia_01.tif"
# var <- "pr"
# month <- "01"

# libraries
library(tidyverse)
library(terra)

# code
base <- rast(base_file)
bias <- rast(bias_file)[[as.numeric(month)]] %>% 
  resample(base)


relative <- FALSE
if(var == "pr")
  relative = TRUE
if(relative)
  warning("Using relative bias.")
if(relative)
  proj <- base + bias
if(!relative)
  proj <- base*(1+bias)
writeRaster(x = proj, filename = fileout)
