import xarray as xr
import rioxarray

geotiff_da = rioxarray.open_rasterio("results/chelsa/raw/pr_01_1979.tif")
geotiff_da.to_netcdf("test.nc")

pr = xr.open_dataset("results/cordex/raw/pr_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_historical/pr_AUS-22_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-HZG-CCLM5-0-15_v1_mon_195001-195012.nc")
tas = xr.open_dataset("results/cordex/raw/tas_AUS-22_MOHC-HadGEM2-ES_CCLM5-0-15_historical/tas_AUS-22_MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-HZG-CCLM5-0-15_v1_mon_195001-195012.nc")

xr.merge([pr, tas])
