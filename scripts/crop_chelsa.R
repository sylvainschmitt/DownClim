# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
filein <- snakemake@input[[1]]
fileout <- snakemake@output[[1]]
foldergadm <- snakemake@input[[2]]
country <- snakemake@params$country

# test
# country <- "New-Caledonia"
# filein <- "results/chelsa/raw/pr_01_1979.tif"
# folderout <- "results/countries/New-Caledonia_gadm/"

# libraries
library(tidyverse)
library(geodata)
library(sf)
library(terra)

# code
country <- gsub("-", " ", country)
bb <- gadm(country, level = 0, path = foldergadm) %>% 
  st_as_sf() %>% 
  st_union() %>% 
  st_bbox() %>% 
  st_as_sfc()
rast(filein) %>% 
  crop(bb) %>% 
  writeRaster(fileout)

