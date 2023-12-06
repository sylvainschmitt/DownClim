# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj =  snakemake.output.proj
base =  snakemake.output.base
period_base = snakemake.params.period_base
period_proj = snakemake.params.period_proj
ds = snakemake.output[0]

# test
# country_file = "results/countries/New-Caledonia.shp"

# libs
import xarray as xr

# code
