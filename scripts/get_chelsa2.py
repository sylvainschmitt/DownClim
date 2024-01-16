# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_files = snakemake.input
areas = snakemake.params.area
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
base_years = snakemake.params.base_years
temp_fold = snakemake.params.tmp
nc_files = snakemake.output
cores = snakemake.threads

# test
# area_files = ["results/countries/New-Caledonia.shp", "results/countries/New-Zealand.shp"]
# areas=["New-Caledonia", "New-Zealand"]
# variables = ["tas", "pr"]
# time_frequency = "mon"
# temp_fold = "test"
# cores = 2
# base_years = "1980-2018"
# nc_files = ["results/chelsa2/raw/New-Caledonia_chelsa2.nc", "results/chelsa2/raw/New-Zealand_chelsa2.nc"]

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
def get_year(year, areas, var, time_freq, month_max, temp_fold):
        mm = month_max + 1
        a = {key: list() for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for month in range(1, mm):
                url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/' + time_freq + '/' + var + "/" \
                        + 'CHELSA_' + var + '_' + '%02d' % (month,) + '_' + str(year) + "_V.2.1.tif"
                ds = rio.open_rasterio(url, decode_coords="all").to_dataset('band').rename_vars({1 : "var"})
                for area in areas: 
                        a[area.NAME_0.values[0]].append(ds.rio.clip_box(minx=area.bounds.minx[0], miny=area.bounds.miny[0],
                                                        maxx=area.bounds.maxx[0], maxy=area.bounds.maxy[0]))
        ds_year = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for area in areas: 
                ds_year[area.NAME_0.values[0]] = xr.concat([i for i in a[area.NAME_0.values[0]]], pd.Index(pd.date_range(datetime.datetime(year, 1, 1), periods=month_max, freq="M"), name="time"))
        del a
        for area in areas: 
                ds_year[area.NAME_0.values[0]] = ds_year[area.NAME_0.values[0]][["time", "x", "y", "var"]]
        paths = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for area in areas: 
                paths[area.NAME_0.values[0]] = temp_fold + "/" + area.NAME_0.values[0] + "_" + var + "_" + str(year) + ".nc"
        for area in areas: 
                ds_year[area.NAME_0.values[0]].to_netcdf(paths[area.NAME_0.values[0]])
        del ds_year
        return(paths)

def get_var(years, areas, var, time_freq, month_max, temp_fold, cores):
        if var == "pr" and years[1] >= 2013:
                years.remove(2013) # issue to heck with the tif
        if var == "pr" and years[1] >= 2016:
                years.remove(2016) # issue to heck with the tif
        pool = mp.Pool(cores)
        paths = pool.starmap_async(get_year, [(y, areas, var, time_freq, month_max, temp_fold) for y in years]).get()
        pool.close()
        paths2 = {key: list() for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for i in list(range(len(years))):
                for area in areas:
                        paths2[area.NAME_0.values[0]].append(paths[i][area.NAME_0.values[0]])
        del paths
        a = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for area in areas:
                a[area.NAME_0.values[0]] = list(map(xr.open_dataset, paths2[area.NAME_0.values[0]]))
        for area in areas:
                for p in paths2[area.NAME_0.values[0]]:
                        os.remove(p)             
        ds_var = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for area in areas:
                ds_var[area.NAME_0.values[0]] = xr.concat(a[area.NAME_0.values[0]], "time")
        del a
        for area in areas:
                if var == "pr":
                        ds_var[area.NAME_0.values[0]]["var"] = ds_var[area.NAME_0.values[0]]["var"] * 0.01
                else:
                        ds_var[area.NAME_0.values[0]]["var"] = ds_var[area.NAME_0.values[0]]["var"] * 0.1
        for area in areas:
                if var == "pr":
                        ds_var[area.NAME_0.values[0]] = ds_var[area.NAME_0.values[0]].rename_vars({"var" : "pr"})
                        ds_var[area.NAME_0.values[0]].pr.attrs = {'standard_name': 'precipitation', 'long_name': 'Monthly precipitation',
                                                                'units': 'mm month-1', 
                                                                'explanation' : 'Precipitation" in the earth\'s atmosphere means precipitation of water in all phases.'}
                if var == "tas":
                        ds_var[area.NAME_0.values[0]] = ds_var[area.NAME_0.values[0]].rename_vars({"var" : "tas"})
                        ds_var[area.NAME_0.values[0]].tas.attrs = {'standard_name': 'temperature at surface', 'long_name': 'Monthly mean daily air temperature',
                                                                'units': 'K', 'explanation' : 'Daily mean air temperatures at 2 metres from hourly ERA5 data.'}
                if var == "tasmin":
                        ds_var[area.NAME_0.values[0]] = ds_var[area.NAME_0.values[0]].rename_vars({"var" : "tasmin"})
                        ds_var[area.NAME_0.values[0]].tasmin.attrs = {'standard_name': 'minimum temperature at surface', 'long_name': 'Monthly minimum daily air temperature',
                                                                'units': 'K', 'explanation' : 'Daily minimum air temperatures at 2 metres from hourly ERA5 data.'}
                if var == "tasmax":
                        ds_var[area.NAME_0.values[0]] = ds_var[area.NAME_0.values[0]].rename_vars({"var" : "tasmax"})
                        ds_var[area.NAME_0.values[0]].tasmax.attrs = {'standard_name': 'maximum temperature at surface', 'long_name': 'Monthly maximum daily air temperature',
                                                                'units': 'K', 'explanation' : 'Daily maximum air temperatures at 2 metres from hourly ERA5 data.'}
        paths = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
        for area in areas: 
                paths[area.NAME_0.values[0]] = temp_fold + "/" + area.NAME_0.values[0] + "_" + var + "_all.nc"
        for area in areas: 
                ds_var[area.NAME_0.values[0]].to_netcdf(paths[area.NAME_0.values[0]])
        del ds_var
        return(paths)

# code
os.mkdir(temp_fold)

areas = list(map(geopandas.read_file, area_files))

base_years = list(map(int, base_years.split("-")))
years = list(range(base_years[0], base_years[1]+1))

if time_frequency == "mon":
        tf = "monthly"
        months = 12
        # months = 2
if time_frequency != "mon":
        raise Exception("Currently only monthly time frequency available!")

paths_var = {key: None for key in variables}
for var in variables:
        paths_var[var] = get_var(years, areas, var, tf, months, temp_fold, cores)        

paths_area = {key: list() for key in list(map(lambda x: x.NAME_0.values[0], areas))}
for area in areas:
        for var in variables:
                paths_area[area.NAME_0.values[0]].append(paths_var[var][area.NAME_0.values[0]])

a = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
for area in areas:
        a[area.NAME_0.values[0]] = list(map(xr.open_dataset, paths_area[area.NAME_0.values[0]]))
ds_all = {key: None for key in list(map(lambda x: x.NAME_0.values[0], areas))}
for area in areas:
        ds_all[area.NAME_0.values[0]] = xr.merge(a[area.NAME_0.values[0]])
del a

for i in list(range(len(nc_files))):
        ds_all[areas[i].NAME_0.values[0]].to_netcdf(nc_files[i])
        
del ds_all
shutil.rmtree(temp_fold)
