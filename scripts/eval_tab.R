# snakemake log
log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file, append = TRUE, type = "message")
sink(log_file, append = TRUE)

# snakemake vars
filein <- snakemake@input[[1]]
fileout <- snakemake@output[[1]]

# libraries
library(tidyverse)

# code
read_tsv(filein) %>% 
  rename(observed = `chelsa raw`) %>% 
  gather(typeorigin, predicted, -variable, -country, -domain, -gcm, -rcm, -rcp, -month, -observed) %>% 
  separate(typeorigin, c("type", "origin")) %>% 
  group_by(variable, country, domain, gcm, rcm, rcp, month, type, origin) %>% 
  summarise(cc = cor(predicted, observed, method = "pearson"), 
            rmse = sqrt(mean((predicted/10-observed/10)^2)),
            bias = mean(predicted/10-observed/10)) %>% 
  write_tsv(fileout)