# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
proj_file =  snakemake.input.proj
base_file =  snakemake.input.base_eval
ds_file =  snakemake.input.ds
pts_file =  snakemake.input.pts
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
baseline = snakemake.params.baseline
period_base = snakemake.params.period_base
period_eval = snakemake.params.period_eval
base_eval = snakemake.params.base_eval
eval = snakemake.output[0]

# test
# proj_file =  'results/projection/raw/New-Caledonia_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1.nc'
# base_file = 'results/chelsa2/raw/New-Caledonia_chelsa2.nc'
# ds_file =  'results/projection/downscaled/New-Caledonia_CORDEX_none_AUS-22_CLMcom-HZG_MOHC-HadGEM2-ES_rcp26_r1i1p1_CCLM5-0-15_v1_chelsa2_2006-2019_1980-2005_bc.nc'
# pts_file =  'results/countries/New-Caledonia_pts.shp'
# period_eval="2006-2019"
# area="New-Caledonia"
# project="CORDEX"
# activity="none"
# domain="AUS-22"
# institute="CLMcom-HZG"
# model="MOHC-HadGEM2-ES"
# experiment="rcp26"
# ensemble="r1i1p1"
# rcm="CCLM5-0-15"
# downscaling="v1"
# baseline="chelsa2"
# period_eval="2006-2019"
# period_base="1980-2005"
# base_eval="chelsa2"

# libs
import xarray as xr
import geopandas
import pandas as pd

# open and extract
pts = geopandas.read_file(pts_file).centroid
x_indexer = xr.DataArray(pts.x, dims=["point"])
y_indexer = xr.DataArray(pts.y, dims=["point"])

# dates 
eval_dmin = period_eval.split("-")[0] + "-01-01" # only january for the moment !
eval_dmax = period_eval.split("-")[1] + "-01-01"

# projection
proj = xr.open_dataset(proj_file).sel(time=slice(eval_dmin, eval_dmax)).groupby("time.month").mean("time").sel(lon=x_indexer, lat=y_indexer, method="nearest")
proj_df = proj.to_dataframe()
proj_df.reset_index(inplace=True)
proj_df = proj_df.drop(['point', "height"], axis = 1).melt(id_vars=["month", "lon", "lat"])
proj_df["type"] = "projection"
proj_df["source"] = project
proj_df["status"] = "raw"
proj_df["method"] = ""
proj_df["area"] = area
proj_df["project"] = project
proj_df["activity"] = activity
proj_df["domain"] = domain
proj_df["institute"] = institute
proj_df["model"] = model
proj_df["experiment"] = experiment
proj_df["ensemble"] = ensemble
proj_df["rcm"] = rcm
proj_df["downscaling"] = downscaling
proj_df["baseline"] = baseline
proj_df["period_eval"] = period_eval
proj_df["period_base"] = period_base
proj_df["base_eval"] = base_eval

# base
base = xr.open_dataset(base_file).sel(time=slice(eval_dmin, eval_dmax)).groupby("time.month").mean("time").rename({"x" : "lon", "y" : "lat"}).sel(lon=x_indexer, lat=y_indexer, method="nearest")
base_df = base.to_dataframe()
base_df.reset_index(inplace=True)
base_df = base_df.drop(['point', "spatial_ref"], axis = 1).melt(id_vars=["month", "lon", "lat"])
base_df["type"] = "baseline"
base_df["source"] = baseline
base_df["status"] = "raw"
base_df["method"] = ""
base_df["area"] = area
base_df["project"] = ""
base_df["activity"] = ""
base_df["domain"] = ""
base_df["institute"] = ""
base_df["model"] = ""
base_df["experiment"] = ""
base_df["ensemble"] = ""
base_df["rcm"] = ""
base_df["downscaling"] = ""
base_df["baseline"] = ""
base_df["period_eval"] = ""
base_df["period_base"] = ""
base_df["base_eval"] = ""

# ds
ds = xr.open_dataset(ds_file).rename({"x" : "lon", "y" : "lat"}).sel(lon=x_indexer, lat=y_indexer, method="nearest")
ds_df = ds.to_dataframe()
ds_df.reset_index(inplace=True)
ds_df = ds_df.drop(['point', "height"], axis = 1).melt(id_vars=["month", "lon", "lat"])
ds_df["type"] = "projection"
ds_df["source"] = project
ds_df["status"] = "downscaled"
ds_df["method"] = "bias correction"
ds_df["area"] = area
ds_df["project"] = project
ds_df["activity"] = activity
ds_df["domain"] = domain
ds_df["institute"] = institute
ds_df["model"] = model
ds_df["experiment"] = experiment
ds_df["ensemble"] = ensemble
ds_df["rcm"] = rcm
ds_df["downscaling"] = downscaling
ds_df["baseline"] = baseline
ds_df["period_eval"] = period_eval
ds_df["period_base"] = period_base
ds_df["base_eval"] = base_eval

# gather & write
pd.concat([proj_df, base_df, ds_df]).to_csv(eval, sep="\t")
