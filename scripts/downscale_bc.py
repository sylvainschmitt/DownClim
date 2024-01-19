# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj_proj_file =  snakemake.input.proj_proj
proj_base_file =  snakemake.input.proj_base
base_base_file =  snakemake.input.base_base
baseline = snakemake.params.base
ds = snakemake.output[0]

# test
# proj_proj_file="results/projection/means/Vanuatu_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1_chelsa2_monthly-means_2006-2019.nc"
# proj_base_file="results/projection/means/Vanuatu_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1_chelsa2_monthly-means_1980-2005.nc"
# base_base_file="results/chelsa2/means/Vanuatu_chelsa2_monthly-means_1980-2005.nc"
# baseline="chelsa2"

# libs
import xarray as xr
import xesmf as xe

# open
proj_proj = xr.open_mfdataset(proj_proj_file, parallel=True)
proj_base = xr.open_mfdataset(proj_base_file, parallel=True)
base_base = xr.open_mfdataset(base_base_file, parallel=True)

# anomalies
anomalies = proj_proj - proj_base
if 'pr' in list(proj_proj.keys()):
  anomalies_rel =  (proj_proj - proj_base)/proj_base
  anomalies["pr"] = anomalies_rel["pr"]
  
# interpolate
regridder = xe.Regridder(anomalies, base_base, "bilinear")
anomalies_ds = regridder(anomalies, keep_attrs=True)

# add to the baseline
proj_ds = base_base + anomalies_ds
if 'pr' in list(proj_proj.keys()):
  proj_ds = base_base * (1 + anomalies_ds)
  proj_ds["pr"] = proj_ds["pr"]
proj_ds = proj_ds.assign_attrs(proj_proj.attrs | anomalies_ds.attrs)
proj_ds.attrs['downscaling'] = "Downscaled with DowClim v0.1.0"
proj_ds.attrs['downscaling_method'] = "Bias correction"
proj_ds.attrs['downscaling_baseline'] = baseline

# write
delayed = proj_ds.to_netcdf(ds, compute=False)
results = delayed.compute(scheduler='threads') # we may need to limit dask to a defined number of cores
