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

# assembly test
# import xarray as xr
# import rioxarray
# 
# geotiff_da = rioxarray.open_rasterio("results/chelsa/raw/pr_01_1979.tif")
# geotiff_da.to_netcdf("test.nc")
# 
# pr = xr.open_dataset("results/cordex/raw/pr_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_historical/pr_AUS-22_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-HZG-CCLM5-0-15_v1_mon_195001-195012.nc")
# tas = xr.open_dataset("results/cordex/raw/tas_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_historical/tas_AUS-22_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-HZG-CCLM5-0-15_v1_mon_195001-195012.nc")
# 
# xr.merge([pr, tas])
