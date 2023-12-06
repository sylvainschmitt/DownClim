# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj_files =  snakemake.input
period_base = snakemake.params.period_base
period_proj = snakemake.params.period_proj
ensemble = snakemake.output[0]

# test
raise Exception("Currently in development!")

# libs
import xarray as xr

# code
