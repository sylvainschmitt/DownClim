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
esgf_credential = snakemake.params.esgf_credential
folder = snakemake.output[0]
cores = snakemake.threads

# test
base_files = ["results/chelsa2/raw/New-Caledonia_chelsa2.nc", "results/chelsa2/raw/Vanuatu_chelsa2.nc"]
folder = "results/projection/raw/CMIP5_world_NCC_NorESM1-M_rcp26_r1i1p1_none_none_chelsa2"
areas = ["New-Caledonia", "Vanuatu"]
institute = "MOHC"
model = "HadGEM2-ES"
experiment = "rcp26"
ensemble = "r1i1p1"
variables = ["tas", "pr"]
time_frequency = "mon"
esgf_credential = "config/credentials_esgf.yml"
cores = 10

# libs
import os
from pyesgf.search import SearchConnection
import yaml
from pyesgf.logon import LogonManager
import xarray as xr
import xesmf as xe
from datetime import datetime as dt
import numpy as np

# funs
def convert_cf_to_dt(x):
  return dt.strptime(str(x), '%Y-%m-%d %H:%M:%S')

# list
server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)
ctx = conn.new_context(
    project="CMIP5",
    institute = institute,
    model = model,
    experiment=[experiment, "historical"],    
    ensemble = ensemble,
    variable=variables,
    time_frequency=time_frequency,
    realm='atmos'
    data_node="esgf.ceda.ac.uk"
    )
# ctx.hit_count
# list(map(lambda res: res.dataset_id, ctx.search()))
all_results = list(map(lambda res: res.file_context().search(), ctx.search()))
all_files = []
for res in all_results: 
    for file in res:
        if any(var in file.opendap_url for var in list(map(lambda x: x + "_", variables))): # the variables filter above is useless!
            all_files.append(file.opendap_url)

# connect
creds = yaml.safe_load(open(esgf_credential, 'r'))
lm = LogonManager()
lm.logon_with_openid(openid=creds['openid'], password=creds['pwd'], interactive=False)

# read & prepare
ds = xr.open_mfdataset(all_files, combine='nested', concat_dim='time', parallel=True)
cf = type(ds["time"].values[0]) is not np.datetime64
if cf: # only cftime if not dt but should include more cases
  ds['time'] = [*map(convert_cf_to_dt, ds.time.values)] 
if 'pr' in list(ds.keys()):
  ds['pr'] = ds.pr * 60*60*24*ds.time.dt.days_in_month  # s-1 to month-1 with 30 days per month
  ds.pr.attrs["units"] = 'mm month-1'

# regrid and write per country
os.mkdir(folder)
for i in list(range(len(areas))):
  base = xr.open_dataset(base_files[i])
  regridder = xe.Regridder(ds, base, "bilinear")
  ds_r = regridder(ds, keep_attrs=True)
  path = folder + "/" + areas[i] + '.nc'
  delayed = ds_r.to_netcdf(path, compute=False)
  results = delayed.compute(scheduler='threads') # we may need to limit dask to a defined number of cores
  