import fsspec
import xarray as xr
import pandas as pd

xmin = 155.86890
xmax = 172.09009
ymin = -22.84806
ymax = -17.39917
variable = "pr"

a = []
for month in range(1, 13):
  url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/ncdf/CHELSA_' + \
        variable + '_' + '%02d' % (month,) + '_1981-2010_V.2.1.nc'
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
ds.to_netcdf("test.nc")

# ref = xr.open_dataset("results/cordex/cropped/pr_New-Caledonia_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_rcp26.nc")
