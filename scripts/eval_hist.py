# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
in_file =  snakemake.input[0]
out_file = snakemake.output[0]

# test
# in_file = "results/projection/downscaled/New-Caledonia_CORDEX_AUS-22_ICTP_NCC-NorESM1-M_rcp26_r1i1p1_RegCM4-7_v0_chelsa2_monthly-means_2006-2019_1980-2005_bc.nc"

# code
import pandas as pd     
import numpy as np   
import xarray as xr
ds = xr.open_dataset(in_file)
# need to mask land, maybe as an option

# define bins per variable with configurable pars    
# to do per month    
bins= np.arange(0, 50, 1)
labels= np.arange(0.5, 49.5, 1)
out = pd.cut(ds.tas.values.ravel(), bins=bins, labels=labels)
res = out.value_counts()
# add a variable name column
# add a model column, need to add variable from wildcard
# bind rows
# write file
# res.to_csv("test.csv")
            