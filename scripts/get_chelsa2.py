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
# area_files = ["results/countries/New-Caledonia.shp", "results/countries/Vanuatu.shp"]
# areas=["New-Caledonia", "Vanuatu"]
# variables = ["tas", "pr"]
# time_frequency = "mon"
# temp_fold = "test"
# cores = 2
# base_years = "1997-1998"
# nc_files = ["results/chelsa2/raw/New-Caledonia_chelsa2.nc", "results/chelsa2/raw/Vanuatu_chelsa2.nc"]

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
        areas_names = list(map(lambda x: x.NAME_0.values[0], areas))
        print("Getting year " + str(year) + "..." + " for areas " + str(areas_names))
        areas_bounds = list(map(lambda x: x.bounds, areas))
        a = {key: list() for key in areas_names}
        for month in range(1, month_max+1):
                url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/' + time_freq + '/' + var + "/" \
                        + 'CHELSA_' + var + '_' + '%02d' % (month,) + '_' + str(year) + "_V.2.1.tif"
                ds = rio.open_rasterio(url, decode_coords="all").to_dataset('band').rename_vars({1 : var})
                for i, _ in enumerate(areas): 
                        a[areas_names[i]].append(ds.rio.clip_box(minx=areas_bounds[i].minx[0], 
                                                                 miny=areas_bounds[i].miny[0],
                                                                 maxx=areas_bounds[i].maxx[0], 
                                                                 maxy=areas_bounds[i].maxy[0]))
        paths = {}
        for area_name in areas_names: 
                ds_year = xr.concat(a[area_name], pd.Index(pd.date_range(datetime.datetime(year, 1, 1), periods=month_max, freq="M"), name="time"))
                del a[area_name]
                ds_year = ds_year[["time", "x", "y", var]]
                paths[area_name] = temp_fold + "/" + area_name + "_" + var + "_" + str(year) + ".nc"
#                ds_year.chunk({'time':1, 'x':400, 'y':400}).to_netcdf(paths[area_name], unlimited_dims="time")
                ds_year.to_netcdf(paths[area_name], unlimited_dims="time")
                del ds_year
        return(paths)

def get_var(years, areas, var, time_freq, month_max, temp_fold, cores):
        print("Getting variable " + var + "...")
        if var == "pr" and years[1] >= 2013:
                years.remove(2013) # issue to heck with the tif
        if var == "pr" and years[1] >= 2016:
                years.remove(2016) # issue to heck with the tif
        pool = mp.Pool(cores)
        paths = pool.starmap_async(get_year, [(y, areas, var, time_freq, month_max, temp_fold) for y in years]).get()
        pool.close()
        areas_names = list(map(lambda x: x.NAME_0.values[0], areas))
        paths2 = {key: list() for key in areas_names}
        paths_concat = {}
        print("Concatenating all years")
        for area_name in areas_names:
                [paths2[area_name].append(path[area_name]) for path in paths]
#                ds_var = xr.concat(ds_var_list[area_name], dim="time")
                ds_var = xr.open_mfdataset(paths2[area_name])
                if var == "pr":
                        ds_var[var].values = ds_var[var].values * 0.01
                        ds_var.pr.attrs = {'standard_name': 'precipitation', 
                                           'long_name': 'Monthly precipitation',
                                           'units': 'mm month-1', 
                                           'explanation' : 'Precipitation" in the earth\'s atmosphere means precipitation of water in all phases.'}
                elif var == "tas":
                        ds_var[var] = ds_var[var] * 0.1
                        ds_var.tas.attrs = {'standard_name': 'temperature at surface', 
                                            'long_name': 'Monthly mean daily air temperature',
                                            'units': 'K', 
                                            'explanation' : 'Daily mean air temperatures at 2 meters.'}
                elif var == "tasmin":
                        ds_var[var] = ds_var[var] * 0.1
                        ds_var.tasmin.attrs = {'standard_name': 'minimum temperature at surface', 
                                               'long_name': 'Monthly minimum daily air temperature',
                                               'units': 'K', 
                                               'explanation' : 'Daily minimum air temperatures at 2 meters.'}
                elif var == "tasmax":
                        ds_var[var] = ds_var[var] * 0.1
                        ds_var.tasmax.attrs = {'standard_name': 'maximum temperature at surface', 
                                               'long_name': 'Monthly maximum daily air temperature',
                                               'units': 'K', 
                                               'explanation' : 'Daily maximum air temperatures at 2 meters.'}
                else:   
                        "Problem, variable "+ var + " not recognized, you should use one of the following: tas, tasmin, tasmax, pr."
                paths_concat[area_name] = temp_fold + "/" + area_name + "_" + var + "_all.nc"
                ds_var.to_netcdf(paths_concat[area_name], encoding = {var: {"dtype": "float32", "zlib": True}})
                del ds_var
                for p in paths2[area_name]:
                        os.remove(p)             
        return(paths_concat)

# code
os.mkdir(temp_fold)

areas = list(map(geopandas.read_file, area_files))
base_years = list(map(int, base_years.split("-")))
years = list(range(base_years[0], base_years[1]+1))
areas_names = list(map(lambda x: x.NAME_0.values[0], areas))

if time_frequency == "mon":
        tf = "monthly"
        months = 12
if time_frequency != "mon":
        raise Exception("Currently only monthly time frequency available!")

paths_var = {}
for var in variables:
        paths_var[var] = get_var(years, areas, var, tf, months, temp_fold, cores)        

paths_areas = {}
for i, area_name in enumerate(areas_names):
        print("Merging files for area " + area_name + "...")
        paths_areas[area_name] = [paths_var[var][area_name] for var in variables]
        xr.merge(list(map(lambda f: xr.open_dataset(f),paths_areas[area_name]))).to_netcdf(nc_files[i])
        for p in paths_areas[area_name]:
                os.remove(p)
shutil.rmtree(temp_fold)
