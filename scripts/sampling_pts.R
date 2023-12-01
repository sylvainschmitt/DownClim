# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
foldergadm <- snakemake@input[[1]]
fileout <- snakemake@output[[1]]
country <- snakemake@params$country
log10_eval_pts <- snakemake@params$log10_eval_pts

# test
# foldergadm <- "results/countries/New-Caledonia_gadm"
# country <- "New-Caledonia"
# log10_eval_pts <- 4

# libraries
library(tidyverse)
library(geodata)
library(terra)
library(sf)

# code
land <- gadm(gsub("-", " ", country), level = 0, path = foldergadm) %>% 
  st_as_sf() %>% 
  st_union()
pts <- st_sample(land, 10^log10_eval_pts) %>% 
  st_coordinates() %>% 
  as.data.frame()
write_tsv(pts, fileout)
