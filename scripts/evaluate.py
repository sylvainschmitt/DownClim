# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj_file =  snakemake.input.proj
base_file =  snakemake.input.base
ds_file =  snakemake.input.ds
shp_file =  snakemake.input.shp
eval = snakemake.output[0]

# test
# proj_file =  'results/projection/means/New-Caledonia_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1_chelsa2_monthly-means_2006-2019.nc'
# base_file = 'results/chelsa2/means/New-Caledonia_chelsa2_monthly-means_2006-2019.nc'
# ds_file =  'results/projection/downscaled/New-Caledonia_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1_chelsa2_monthly-means_2006-2019_1980-2005_bc.nc'
# shp_file =  'results/countries/New-Caledonia.shp'

import rioxarray as rio
import geopandas as gp
import xarray as xr
geodf = gp.read_file(shp_file)
base = rio.open_rasterio(base_file).rio.write_crs("epsg:4326", inplace=True).rio.clip(geodf.geometry.values).to_dataset().expand_dims(dim={"type": ["raw"]}, axis=0)
proj = rio.open_rasterio(proj_file).rio.write_crs("epsg:4326", inplace=True).rio.clip(geodf.geometry.values).to_dataset().expand_dims(dim={"type": ["projection"]}, axis=0)
ds = rio.open_rasterio(ds_file).rio.write_crs("epsg:4326", inplace=True).rio.clip(geodf.geometry.values).to_dataset().expand_dims(dim={"type": ["downscaled"]}, axis=0)
total = ds.combine_first(proj).combine_first(base).assign_attrs(xr.open_dataset(ds_file).attrs)
total = total.assign_attrs(xr.open_dataset(ds_file).attrs)
total.to_netcdf(eval)
