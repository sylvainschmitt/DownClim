# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.input[0]
domain = snakemake.params.domain
area = snakemake.params.area
project = snakemake.params.project
activity = snakemake.params.activity
domain = snakemake.params.domain
institute = snakemake.params.institute
model = snakemake.params.model
experiment = snakemake.params.experiment
ensemble = snakemake.params.ensemble
rcm = snakemake.params.rcm
downscaling = snakemake.params.downscaling
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
esgf_credential = snakemake.params.esgf_credential
nc_file = snakemake.output[0]
cores = snakemake.threads

# test
# area_file = "results/countries/French-Guiana.shp"
# nc_file = "test.nc"
# area = "French-Guiana"
# project = "CORDEX"
# activity = "none"
# domain = "AFR-22"
# institute = "CLMcom-HZG"
# model = "NCC-NorESM1-M"
# experiment = "rcp85"
# ensemble = "r1i1p1"
# rcm = "CCLM-0-15"
# downscaling = "v1"
# variables = ["tas", "tasmin", "tasmax", "pr"]
# time_frequency = "mon"
# proj_years = "2071-2100"
# esgf_credential = "config/credentials_esgf.yml"
# cores = 10

# libs
from pyesgf.search import SearchConnection
import yaml
from pyesgf.logon import LogonManager
import xarray as xr
import geopandas
import nctoolkit as nc
import dask.multiprocessing
import re
from datetime import datetime as dt
import numpy

# funs
def convert_to_dt(x):
    return dt.strptime(str(x), '%Y-%m-%d %H:%M:%S')

# list
server = 'https://esgf-node.ipsl.upmc.fr/esg-search/'
conn = SearchConnection(server, distrib=True)
# list of keys: https://esgf-node.llnl.gov/esg-search/search?project=CORDEX&facets=*&limit=0
if(project == "CORDEX"):
  ctx = conn.new_context(
    facets='*',
    project = project,
    domain = domain,
    institute = institute,
    driving_model = model,
    experiment = [experiment, "historical"],
    ensemble = ensemble,
    rcm_name= rcm,
    rcm_version = downscaling,
    time_frequency = time_frequency,
    variable = variables
    )  
if(project != "CORDEX"):
  raise Exception("Currently in development!")
results = ctx.search()
datasets = [res.dataset_id for res in results]
all_results = list(map(lambda res: res.file_context().search(), results))
all_files = []
for res in all_results: 
  all_files = all_files + [file.opendap_url for file in res]

# connect
creds = yaml.safe_load(open(esgf_credential, 'r'))
lm = LogonManager()
lm.logon_with_openid(openid=creds['openid'], password=creds['pwd'], interactive=False)
lm.is_logged_on()

# read, prepare, and write
area = geopandas.read_file(area_file)
res = int(re.findall(r'\d+', domain)[0])/100
ds = xr.open_mfdataset(all_files, parallel=True)
if type(ds["time"].values[0]) is not numpy.datetime64:
  ds['time'] = [*map(convert_to_dt, ds.time.values)]
if 'pr' in list(ds.keys()):
  ds['pr'] = ds['pr']*60*60*24*30 # s-1 to month-1
  ds.pr.attrs["units"] = 'mm month-1'
ds = ds.compute(scheduler='threads')
ds_nc = nc.from_xarray(ds)
ds_nc.to_latlon(lon = [area.bounds.minx[0], area.bounds.maxx[0]], lat = [area.bounds.miny[0], area.bounds.maxy[0]], res = res)
ds = ds_nc.to_xarray()
ds.to_netcdf(nc_file)

# from matplotlib import pyplot as plt # check with plot
# ds.pr[0].plot()
# plt.show()
