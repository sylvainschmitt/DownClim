# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
ds =  snakemake.input.ds
raw =  snakemake.input.raw
base =  snakemake.input.base
out_file = snakemake.output[0]

# test
# ds = ["results/evaluation/hist/ds/New-Caledonia_CORDEX_AUS-22_ICTP_NCC-NorESM1-M_rcp26_r1i1p1_RegCM4-7_v0_chelsa2_monthly-means_2006-2019_1980-2005_bc.tsv", 
#       "results/evaluation/hist/ds/Vanuatu_CORDEX_AUS-22_ICTP_NCC-NorESM1-M_rcp26_r1i1p1_RegCM4-7_v0_chelsa2_monthly-means_2006-2019_1980-2005_bc.tsv", 
#       "results/evaluation/hist/ds/New-Caledonia_CMIP6_world_MRI_MRI-ESM2-0_ssp126_r1i1p1f1_none_none_chelsa2_monthly-means_2006-2019_1980-2005_bc.tsv", 
#       "results/evaluation/hist/ds/Vanuatu_CMIP6_world_MRI_MRI-ESM2-0_ssp126_r1i1p1f1_none_none_chelsa2_monthly-means_2006-2019_1980-2005_bc.tsv"]
# raw = ["results/evaluation/hist/raw/New-Caledonia_CORDEX_AUS-22_ICTP_NCC-NorESM1-M_rcp26_r1i1p1_RegCM4-7_v0_chelsa2_monthly-means_2006-2019.tsv", 
#        "results/evaluation/hist/raw/Vanuatu_CORDEX_AUS-22_ICTP_NCC-NorESM1-M_rcp26_r1i1p1_RegCM4-7_v0_chelsa2_monthly-means_2006-2019.tsv", 
#        "results/evaluation/hist/raw/New-Caledonia_CMIP6_world_MRI_MRI-ESM2-0_ssp126_r1i1p1f1_none_none_chelsa2_monthly-means_2006-2019.tsv", 
#        "results/evaluation/hist/raw/Vanuatu_CMIP6_world_MRI_MRI-ESM2-0_ssp126_r1i1p1f1_none_none_chelsa2_monthly-means_2006-2019.tsv"]
# base = ["results/evaluation/hist/base/New-Caledonia_chelsa2_monthly-means_2006-2019.tsv", 
#         "results/evaluation/hist/base/Vanuatu_chelsa2_monthly-means_2006-2019.tsv"]

# libs
import pandas as pd 

# code
ds_all = pd.concat(list(map(lambda t: pd.read_csv(t, sep="\t"), ds)))
raw_all = pd.concat(list(map(lambda t: pd.read_csv(t, sep="\t"), raw)))
base_all = pd.concat(list(map(lambda t: pd.read_csv(t, sep="\t"), base)))
all = pd.concat([ds_all, raw_all, base_all])
all = all[all["count"] > 0]
all.to_csv(out_file, sep="\t", index=False)
