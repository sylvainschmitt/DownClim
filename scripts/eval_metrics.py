# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
in_file =  snakemake.input[0]
base_file = snakemake.input[1]
area_file =  snakemake.input[2]
out_file = snakemake.output[0]
area = snakemake.params.area
origin = snakemake.params.origin
type = snakemake.params.type
domain = snakemake.params.domain
institute = snakemake.params.institute
model = snakemake.params.model
experiment = snakemake.params.experiment
ensemble = snakemake.params.ensemble
rcm = snakemake.params.rcm
downscaling = snakemake.params.downscaling
base = snakemake.params.base
aggregation = snakemake.params.aggregation
period_proj = snakemake.params.period_proj
period_proj = snakemake.params.period_proj
period_eval = snakemake.params.period_eval
ds_method = snakemake.params.ds_method
base_eval = snakemake.params.base_eval
      
# test
# in_file = "results/projection/means/New-Caledonia_CMIP6_world_MRI_MRI-ESM2-0_ssp126_r1i1p1f1_none_none_chelsa2_monthly-means_2006-2019.nc"
# base_file = "results/chelsa2/means/New-Caledonia_chelsa2_monthly-means_2006-2019.nc"
# area_file = "results/countries/New-Caledonia.shp"

# libs
import pandas as pd     
import numpy as np   
import numpy.ma as ma
import xarray as xr
import geopandas as gp

# code
area_shp = gp.read_file(area_file)
ds = xr.open_dataset(in_file).rio.write_crs("epsg:4326", inplace=True).rio.clip(area_shp.geometry.values, area_shp.crs) # crs should be checked earlier, concern only CORDEX raw
base = xr.open_dataset(base_file).rio.clip(area_shp.geometry.values, area_shp.crs)
variables = list(ds.keys())
months = list(range(1,13))
a = []
for v in variables:
    for m in months:
        pred = ds.sel(month=m)[v].values.ravel()
        pred = pred[~np.isnan(pred)]
        obs = base.sel(month=m)[v].values.ravel()
        obs = obs[~np.isnan(obs)]
        d = {
            'metric': ['CC', 'RMSE', "SDE", "bias"],
            'value': [np.corrcoef(pred, obs)[1,0], np.sqrt(np.mean(pow(pred - obs, 2))), np.std(pred - obs), np.mean(pred - obs)]
            }
        res = pd.DataFrame(data = d)
        res.insert(0, "month", m)
        res.insert(0, "variable", v)
        a.append(res)
tab = pd.concat(a)
tab.insert(0, "area", area)
tab.insert(0, "origin", origin)
tab.insert(0, "type", type)
tab.insert(0, "domain", domain)
tab.insert(0, "institute", institute)
tab.insert(0, "model", model)
tab.insert(0, "experiment", experiment)
tab.insert(0, "ensemble", ensemble)
tab.insert(0, "rcm", rcm)
tab.insert(0, "downscaling", downscaling)
tab.insert(0, "base", base)
tab.insert(0, "aggregation", aggregation)
tab.insert(0, "period_proj", period_proj)
tab.insert(0, "period_eval", period_eval)
tab.insert(0, "ds_method", ds_method)
tab.insert(0, "base_eval", base_eval)
tab[["area", "origin", "type", "domain", "institute", "model", "experiment", "ensemble", "rcm", "downscaling", "base", 
     "aggregation", "period_proj", "period_eval", "ds_method", "base_eval",
     "month", "variable", "metric", "value"]].to_csv(out_file, sep="\t", index=False)
