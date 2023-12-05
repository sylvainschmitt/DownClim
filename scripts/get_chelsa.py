# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
country_file = snakemake.intput[0]
country = snakemake.params.country
var = snakemake.params.var
nc_file = snakemake.output[0]

# test
# country_file = "results/countries/New-Caledonia.shp"
# var = "pr"

# libs
import fsspec
import xarray as xr
import pandas as pd
import geopandas

# code
lim = geopandas.read_file(country_file)
xmin = lim.bounds.minx
xmax = lim.bounds.maxx
ymin = lim.bounds.miny
ymax = lim.bounds.maxx
a = []
for month in range(1, 13):
  url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/ncdf/CHELSA_' + \
        var + '_' + '%02d' % (month,) + '_1981-2010_V.2.1.nc'
  fobj = fsspec.open(url)
  with fsspec.open(url) as fobj:
    ds = xr.open_dataset(fobj).chunk({'lat': 500, 'lon': 500})
    mask_lon = (ds.lon >= xmin) & (ds.lon <= xmax)
    mask_lat = (ds.lat >= ymin) & (ds.lat <= ymax)
    ds = ds.where(mask_lon & mask_lat, drop = True)
    ds.load()
  a.append(ds)
ds = xr.concat([i for i in a], pd.Index([*range(1, 13)], name="month"))
ds = ds.rename({'Band1': 'pr'})
ds.pr.attrs = {'standard_name': 'precipitation', 'long_name': 'Monthly precipitation', 'units': 'kg m-2 month-1 / 10'}
ds.month.attrs = {'standard_name': 'month', 'long_name': 'month'}
ds.to_netcdf(nc_file)

