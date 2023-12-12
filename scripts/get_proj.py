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
proj_years = snakemake.params.proj_years
esgf_credential = snakemake.params.esgf_credential
nc_file = snakemake.output[0]
cores = snakemake.threads

# test
# area_file = "results/countries/New-Caledonia.shp"
# nc_file = "test.nc"
# area = "New-Caledonia"
# project = "CORDEX"
# activity = "none"
# domain = "AUS-22"
# institute = "CLMcom-HZG"
# model = "MOHC-HadGEM2-ES"
# experiment = "rcp26"
# ensemble = "r1i1p1"
# rcm = "CCLM5-0-15"
# downscaling = "v1"
# variables = ["tas", "tasmin", "tasmax", "pr"]
# time_frequency = "mon"
# proj_years = "2071-2100"
# esgf_credential = "config/credentials_esgf.yml"
# cores = 20

# libs
from pyesgf.search import SearchConnection
import yaml
from pyesgf.logon import LogonManager
import xarray as xr
import geopandas
import nctoolkit as nc
import dask.multiprocessing
import re

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
    experiment = experiment,
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
print(all_files)

# connect
creds = yaml.safe_load(open(esgf_credential, 'r'))
lm = LogonManager()
lm.logon_with_openid(openid=creds['openid'], password=creds['pwd'], interactive=False)
lm.is_logged_on()

# read, prepare, and write
area = geopandas.read_file(area_file)
res = int(re.findall(r'\d+', domain)[0])/100
ds = xr.open_mfdataset(all_files, parallel=True)
ds = ds.compute(scheduler='threads')
ds_nc = nc.from_xarray(ds)
ds_nc.to_latlon(lon = [area.bounds.minx[0], area.bounds.maxx[0]], lat = [area.bounds.miny[0], area.bounds.maxy[0]], res = res)
ds = ds_nc.to_xarray()
ds.to_netcdf(nc_file)

# from matplotlib import pyplot as plt # check with plot
# ds.pr[0].plot()
# plt.show()
