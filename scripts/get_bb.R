# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
fileout <- snakemake@output[[1]]
folderout <- snakemake@output[[2]]
figout <- snakemake@output[[3]]
country <- snakemake@params$country

# test
# country <- "New-Caledonia"
# folderout <- "test"

# libraries
library(tidyverse)
library(geodata)
library(sf)

# code
country <- gsub("-", " ", country)
polygon <- gadm(country, level = 0, path = folderout)
g <- ggplot(data = st_as_sf(polygon)) +
  geom_sf() +
  theme_void()
ggsave(plot = g, filename = figout, bg = "white")
bb <- st_bbox(polygon)
data.frame(
  var = paste0(names(bb), "=", as.vector(bb))
) %>% write_tsv(col_names = FALSE, file = fileout)

