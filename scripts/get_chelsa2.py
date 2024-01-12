# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.input[0]
area = snakemake.params.area
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
base_years = snakemake.params.base_years
temp_fold = snakemake.params.tmp
nc_file = snakemake.output[0]
cores = snakemake.threads

# test
# area_file = "results/countries/New-Caledonia.shp"
# temp_fold = "test"
# area="New-Caledonia"
# variables = ["tas", "pr"]
# time_frequency = "mon"
# base_years = "1980-2018"
# nc_file = "test.nc"
# cores = 2

# libs
import os
import shutil
import xarray as xr
import rioxarray as rio
import pandas as pd
import geopandas
import datetime
import multiprocessing as mp

# function
def get_year(year, area, var, time_freq, month_max, temp_fold):
       mm = month_max + 1
       a = []
       for month in range(1, mm):
                url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/' + time_freq + '/' + var + "/" \
                        + 'CHELSA_' + var + '_' + '%02d' % (month,) + '_' + str(year) + "_V.2.1.tif"
                ds = rio.open_rasterio(url, decode_coords="all").to_dataset('band').rename_vars({1 : "var"})
                ds = ds.rio.clip_box(minx=area.bounds.minx[0], miny=area.bounds.miny[0],
                                     maxx=area.bounds.maxx[0], maxy=area.bounds.maxy[0])
                a.append(ds)
       ds_year = xr.concat([i for i in a], pd.Index(pd.date_range(datetime.datetime(year, 1, 1), periods=month_max, freq="M"), name="time"))
       del a
       ds_year = ds_year[["time", "x", "y", "var"]]
       path = temp_fold + "/" + var + "_" + str(year) + ".nc"
       ds_year.to_netcdf(path)
       del ds_year
       return(path)

def get_var(years, area, var, time_freq, month_max, temp_fold, cores):
        if var == "pr" and years[1] <= 2013: 
                years.remove(2013) # issue to heck with the tif
        if var == "pr" and years[1] <= 2016:
                years.remove(2016) # issue to heck with the tif
        pool = mp.Pool(cores)
        paths = pool.starmap_async(get_year, [(y, area, var, time_freq, month_max, temp_fold) for y in years]).get()
        pool.close()
        a = list(map(xr.open_dataset, paths))
        for p in paths:
                os.remove(p)
        ds_var = xr.concat(a, "time")
        del a
        ds_var["var"] = ds_var["var"] * 0.1 # pr need to be divided by 100!
        if var == "pr":
                ds_var = ds_var.rename_vars({"var" : "pr"})
                ds_var.pr.attrs = {'standard_name': 'precipitation', 'long_name': 'Monthly precipitation',
                                   'units': 'mm month-1', 
                                   'explanation' : 'Precipitation" in the earth\'s atmosphere means precipitation of water in all phases.'}
        if var == "tas":
                ds_var = ds_var.rename_vars({"var" : "tas"})
                ds_var.tas.attrs = {'standard_name': 'temperature at surface', 'long_name': 'Monthly mean daily air temperature',
                                    'units': 'K', 'explanation' : 'Daily mean air temperatures at 2 metres from hourly ERA5 data.'}
        if var == "tasmin":
                ds_var = ds_var.rename_vars({"var" : "tasmin"})
                ds_var.tasmin.attrs = {'standard_name': 'minimum temperature at surface', 'long_name': 'Monthly minimum daily air temperature',
                                       'units': 'K', 'explanation' : 'Daily minimum air temperatures at 2 metres from hourly ERA5 data.'}
        if var == "tasmax":
                ds_var = ds_var.rename_vars({"var" : "tasmax"})
                ds_var.tasmax.attrs = {'standard_name': 'maximum temperature at surface', 'long_name': 'Monthly maximum daily air temperature',
                                       'units': 'K', 'explanation' : 'Daily maximum air temperatures at 2 metres from hourly ERA5 data.'}
        path = temp_fold + "/" + var + "_all.nc"
        ds_var.to_netcdf(path)
        del ds_var
        return(path)

# code
os.mkdir(temp_fold)

area = geopandas.read_file(area_file)

base_years = list(map(int, base_years.split("-")))
years = list(range(base_years[0], base_years[1]+1))

if time_frequency == "mon":
        tf = "monthly"
        months = 12
        # months = 2
if time_frequency != "mon":
        raise Exception("Currently only monthly time frequency available!")

paths = []
for var in variables:
        path = get_var(years, area, var, tf, months, temp_fold, cores)
        paths.append(path)
a = list(map(xr.open_dataset, paths))
ds_all = xr.merge(a)
del a
ds_all.to_netcdf(nc_file)
del ds_all
shutil.rmtree(temp_fold)
