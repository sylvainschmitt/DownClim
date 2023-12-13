# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj_file =  snakemake.input.proj
base_file =  snakemake.input.base
period_base = snakemake.params.period_base
period_proj = snakemake.params.period_proj
ds = snakemake.output[0]

# test
# proj_file="results/projection/raw/New-Caledonia_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1.nc"
# base_file="results/chelsa2/raw/New-Caledonia_chelsa2.nc"
# period_proj="2071-2100"
# period_base="1980-2005"

# libs
import xarray as xr
import xesmf as xe

# dates 
ref_dmin = period_base.split("-")[0] + "-01-01"
ref_dmax = period_base.split("-")[1] + "-01-01"
proj_dmin = period_proj.split("-")[0] + "-01-01"
proj_dmax = period_proj.split("-")[1] + "-01-01"

# open
base = xr.open_dataset(base_file)
proj = xr.open_dataset(proj_file)

# slice and summarise
base_ref = base.sel(time=slice(ref_dmin, ref_dmax)).groupby("time.month").mean("time")
proj_ref = proj.sel(time=slice(ref_dmin, ref_dmax)).groupby("time.month").mean("time")
proj_proj = proj.sel(time=slice(proj_dmin, proj_dmax)).groupby("time.month").mean("time")

# anomalies
anomalies = proj_proj - proj_ref
if 'pr' in list(proj.keys()):
  anomalies_rel =  (proj_proj - proj_ref)/proj_ref
  anomalies["pr"] = anomalies_rel["pr"]
  
# interpolate with nctoolkit or other
regridder = xe.Regridder(anomalies, base_ref, "bilinear")
anomalies_ds = regridder(anomalies)
anomalies_ds

# add to the baseline
proj_ds = base_ref + anomalies_ds
if 'pr' in list(proj.keys()):
  proj_ds = base_ref * (1 + anomalies_ds)
  proj_ds["pr"] = proj_ds["pr"]
proj_ds.attrs = proj.attrs
proj_ds.attrs['downscaling'] = "Downscaled with DowClim v0.1.0"
proj_ds.attrs['downscaling_method'] = "Bias correction"
proj_ds.attrs['downscaling_baseline'] = "" # to be added to the baseline attrs

# write
proj_ds.to_netcdf(ds)

# from matplotlib import pyplot as plt # check with plot
# proj_ds.pr[0].plot()
# plt.show()
