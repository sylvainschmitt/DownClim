# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.intput[0]
area = snakemake.params.country
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
base_years = snakemake.params.base_years
nc_file = snakemake.output[0]

# test
raise Exception("Currently in development!")

# libs
import xarray as xr
