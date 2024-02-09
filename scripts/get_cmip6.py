# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
base_files = snakemake.input
areas = snakemake.params.area
institute = snakemake.params.institute
model = snakemake.params.model
experiment = snakemake.params.experiment
ensemble = snakemake.params.ensemble
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
folder = snakemake.output[0]
threads = snakemake.threads

# test
# base_files = ["results/chelsa2/raw/New-Caledonia_chelsa2.nc", "results/chelsa2/raw/Vanuatu_chelsa2.nc"]
# folder = "results/projection/raw/CMIP6_world_NCC_NorESM2-MM_ssp126_r1i1p1f1_none_none_chelsa2"
# areas = ["New-Caledonia", "Vanuatu"]
# institute = "MRI"
# model = "MRI-ESM2-0"
# experiment = "ssp126"
# ensemble = "r1i1p1f1"
# variables =  ["tas"]
# time_frequency = "mon"
# threads = 10

# libs
import os
import pandas as pd
import gcsfs
import xarray as xr
import xesmf as xe
import numpy as np
from datetime import datetime as dt

# funs
def convert_cf_to_dt(x):
  return dt.strptime(str(x), '%Y-%m-%d %H:%M:%S')

# conversions
if time_frequency == "mon":
        table_id = "Amon"
else :
    raise Exception("Currently only monthly time frequency available!")

# code
gcs = gcsfs.GCSFileSystem(token='anon')
df = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')
a = []
for exp in ["historical", experiment]:
  if exp=="historical":
    activity="CMIP"
  else:
    activity="ScenarioMIP"
  for var in variables:
    search_string = "activity_id == '" + activity + "' & table_id == '" + table_id + "' & variable_id == '" + var + "' & experiment_id == '" + exp + "' & institution_id == '" + institute + "' & source_id == '" + model + "' & member_id == '" + ensemble + "'"
    df_ta = df.query(search_string)
    zstore = df_ta.zstore.values[-1]
    mapper = gcs.get_mapper(zstore)
    a.append(xr.open_zarr(mapper, consolidated=True))
ds = xr.merge(a)
ds['time'] = np.sort(ds['time'].values)
ds = ds.chunk(chunks = {'time':100, 'lat': 400, 'lon': 400})
cf = type(ds["time"].values[0]) is not np.datetime64
if cf: # only cftime if not dt but should include more cases
  ds['time'] = [*map(convert_cf_to_dt, ds.time.values)] 
if 'pr' in list(ds.keys()):
  ds['pr'] = ds.pr * 60*60*24*ds.time.dt.days_in_month  # s-1 to month-1
  ds.pr.attrs["units"] = 'mm month-1'

os.mkdir(folder)
for i in list(range(len(areas))):
  base = xr.open_dataset(base_files[i])
  regridder = xe.Regridder(ds, base, "bilinear")
  ds_r = regridder(ds, keep_attrs=True)
  path = folder + "/" + areas[i] + '.nc'
  delayed = ds_r.to_netcdf(path, compute=False)
  results = delayed.compute(scheduler='threads') # we may need to limit dask to a defined number of cores, see in this case "from distributed import Client"
