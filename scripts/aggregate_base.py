# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
base_file =  snakemake.input[0]
agg_file =  snakemake.output[0]
period = snakemake.params.period
aggregation = snakemake.params.aggregation

# test
# base_file="results/chelsa2/raw/New-Caledonia_chelsa2.nc"
# agg_file="test.nc"
# period="2006-2019"
# aggregation="monthly-means"

# libs
import xarray as xr

# aggregation
if(aggregation != "monthly-means"):
  raise Exception("Currently only monthly-means available!")

# code 
dmin = period.split("-")[0] + "-01-01"
dmax = period.split("-")[1] + "-01-01"
base = xr.open_mfdataset(base_file, parallel=True).sel(time=slice(dmin, dmax)).groupby("time.month").mean("time")
delayed = base.to_netcdf(agg_file, compute=False)
results = delayed.compute(scheduler='threads') # we may need to limit dask to a defined number of cores
