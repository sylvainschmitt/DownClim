# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj =  snakemake.output.proj
base =  snakemake.output.base
ds =  snakemake.output.ds
pts =  snakemake.output.pts
period_base = snakemake.params.period_base
period_eval = snakemake.params.period_eval
eval = snakemake.output[0]

# test
# country_file = "results/countries/New-Caledonia.shp"

# libs
import xarray as xr

# code
