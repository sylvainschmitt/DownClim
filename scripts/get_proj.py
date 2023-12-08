# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.intput[0]
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

# test
area_file = "results/countries/New-Caledonia.shp"
nc_file = "test.nc"
area = "New-Caledonia"
project = "CORDEX"
activity = "none"
domain = "AUS-22"
institute = "CLMcom-HZG"
model = "MOHC-HadGEM2-ES"
experiment = "rcp26"
ensemble = "r1i1p1"
rcm = "CCLM5-0-15"
downscaling = "v1"
variables = ["tas", "tasmin", "tasmax", "pr"]
time_frequency = "mon"
proj_years = "2071-2100"
esgf_credential = "config/credentials_esgf.yml"

# code
from pyesgf.search import SearchConnection
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

# ctx.hit_count
results = ctx.search()
datasets = [res.dataset_id for res in results]

all_results = []
for res in results:
  print(res.dataset_id)
  files = res.file_context().search()
  all_results.append(files)

all_files = []
for res in all_results:
  all_files = all_files + [file.opendap_url for file in res]

import yaml
creds = yaml.safe_load(open(esgf_credential, 'r'))

from pyesgf.logon import LogonManager
lm = LogonManager()
lm.logon_with_openid(openid=creds['openid'], password=creds['pwd'], interactive=False)
lm.is_logged_on()

import xarray as xr
ds = xr.open_dataset(all_files[0], decode_cf=True, decode_coords="all") # good until here

# import pyproj
# lonmesh, latmesh = np.meshgrid(ds.rlon, ds.rlat)
# source_crs = pyproj.CRS(ds.rio.crs) # Coordinate system of the file
# target_crs = pyproj.CRS(init="epsg:4326") # Global lat-lon coordinate system
# polar_to_latlon = pyproj.Transformer.from_crs(source_crs, target_crs)
# lon, lat = polar_to_latlon.transform(lonmesh, latmesh)
# lon = xr.DataArray(lon, dims=('rlat','rlon'))
# lat = xr.DataArray(lat, dims=('rlat','rlon'))
# ds_lonlat = ds.interp({'rlon':lon, 'rlat':lat}, method='nearest').load()
# ds_lonlat = ds_lonlat.assign_coords({'lon': (('rlon'), ds_lonlat.rlon.data), 'lat': (('rlat'), ds_lonlat.rlat.data)}).drop_vars(['rlon','rlat']).reset_coords()
# ds2 = ds.assign_coords({'lon': (('lon'), ds), 'lat': (('lat'), ds)}).drop_vars(['rlon','rlat']).reset_coords()

import geopandas 
country = geopandas.read_file(area_file)
clipped = ds.rio.clip(country.geometry) # not clipping because wrong projection
clipped.to_netcdf("test.nc")

# import cordex as cx # test with cordex
# from pyproj import CRS
# ds_reproj = cx.transform_coords(ds, src_crs=CRS.from_cf(ds.cf["grid_mapping"].attrs),  trg_crs=CRS.from_string("EPSG:4326")) # not clipping because wrong projection
# clipped = ds_reproj.rio.clip(country.geometry)

# import rioxarray # test with roxarray
# ds.rio.reproject("EPSG:4326")
# ds.rio.estimate_utm_crs()
# ds_utm = ds.rio.reproject(ds.rio.estimate_utm_crs())
# ds_lonlat = ds.rio.reproject("EPSG:4326")
# ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
# ds.rio.write_crs("epsg:4326", inplace=True)

# import xesmf as xe  # test with xesmf
# regridder = xe.Regridder(ds, xe.util.grid_global(0.22, 0.22), "bilinear", periodic=True)
# regridder
# remap_xe = regridder(ds)
# remap_xe.pr[0].plot()

import xarray as xr
ds = xr.open_dataset(all_files[0], decode_cf=True, decode_coords="all")
ds_lonlat = ds.interp(rlon=ds.lon, rlat=ds.lat)
ds_lonlat = ds_lonlat.drop_vars(["lat", "lon"])
ds_lonlat = ds_lonlat.rename_vars({"rlat": "lat", "rlon": "lon"})
ds_lonlat.to_netcdf("test.nc")

from matplotlib import pyplot as plt # check with plot
ds_lonlat.pr[0].plot()
plt.show()
