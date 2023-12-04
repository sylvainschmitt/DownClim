import xarray as xr
import rioxarray

geotiff_da = rioxarray.open_rasterio("results/chelsa/raw/pr_01_1979.tif")
geotiff_da.to_netcdf("test.nc")
