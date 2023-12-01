# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
ds_file <- snakemake@input[[1]]
raw_file <- snakemake@input[[2]]
base_file <- snakemake@input[[3]]
pts_file <- snakemake@input[[4]]
fileout <- snakemake@output[[1]]
var <- snakemake@params$var
month <- snakemake@params$month
country <- snakemake@params$country

# test
# ds_file <- "results/downscaled/eval_pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26_01.tif"
# raw_file <- "results/cordex/eval/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc"
# base_file <- "results/chelsa/eval/pr_New-Caledonia_01.tif"
# pts_file <- "results/countries/New-Caledonia_sampling.tsv"
# var <- "pr"
# month <- "01"
# country <- "New-Caledonia"

# libraries
library(tidyverse)
library(terra)
library(sf)

# infos
# chelsa units
# pr: kg m-2 month-1 / 100
# tas: K / 10
# cordex unites
# pr: kg m-2 s-1
# tas: K

# code
pts <- read_tsv(pts_file)
base <- rast(base_file)
if(var == "pr")
  base <- base/100/(31*24*60*60) # depend on the month!!
if(var != "pr")
  base <- base/10
ds <- rast(ds_file)
if(var == "pr")
  ds <- ds/100/(31*24*60*60) # depend on the month!!
if(var != "pr")
  ds <- ds/10
raw <- rast(raw_file)[[as.numeric(month)]]

infos <- str_split_1(ds_file, "_")
res <- tibble(
  variable = infos[2],
  country = infos[3],
  domain = infos[4],
  gcm = infos[5],
  rcm = infos[6],
  rcp = infos[7],
  month = month,
  "chelsa raw" = extract(base, pts)[,2],
  "cordex raw" = extract(raw, pts)[,2],
  "cordex downscaled" = extract(ds, pts)[,2]
)
write_tsv(res, fileout)
