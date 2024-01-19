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

# test
# area_files = ["results/countries/New-Caledonia.shp", "results/countries/Vanuatu.shp"]
# areas=["New-Caledonia", "Vanuatu"]
# variables = ["tas", "tasmin"]
# time_frequency = "mon"
# temp_fold = "results/chelsa2/raw/tmp"
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

# funs
def get_year(year, areas, var, time_freq, month_max, temp_fold):
        areas_names = list(map(lambda x: x.NAME_0.values[0], areas))
        print("Getting year " + str(year) + " for variables " + var + " and areas " + str(areas_names))
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
                ds_year.chunk({'time':1, 'x':100, 'y':100}).to_netcdf(paths[area_name])
                del ds_year
        return(paths)

def prep_netcdf(ds):
        if 'pr' in list(ds.keys()):
                ds['pr'].values = ds['pr'].values * 0.01
                ds.pr.attrs = {'standard_name': 'precipitation', 
                        'long_name': 'Monthly precipitation',
                        'units': 'mm month-1', 
                        'explanation' : 'Precipitation" in the earth\'s atmosphere means precipitation of water in all phases.'}
        if 'tas' in list(ds.keys()):
                ds['tas'].values = ds['tas'].values * 0.1
                ds.tas.attrs = {'standard_name': 'temperature at surface', 
                        'long_name': 'Monthly mean daily air temperature',
                        'units': 'K', 
                        'explanation' : 'Daily mean air temperatures at 2 meters.'}        
        if 'tasmin' in list(ds.keys()):
                ds['tasmin'].values = ds['tasmin'].values * 0.1
                ds.tasmin.attrs = {'standard_name': 'minimum temperature at surface', 
                        'long_name': 'Monthly minimum daily air temperature',
                        'units': 'K', 
                        'explanation' : 'Daily minimum air temperatures at 2 meters.'}   
        if 'tasmax' in list(ds.keys()):
                ds['tasmax'].values = ds['tasmax'].values * 0.1
                ds.tasmax.attrs = {'standard_name': 'maximum temperature at surface', 
                        'long_name': 'Monthly maximum daily air temperature',
                        'units': 'K', 
                        'explanation' : 'Daily maximum air temperatures at 2 meters.'}
        return(ds)
  
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

paths = {key: list() for key in areas_names}
for v in variables:
        for y in years:
                raw_paths = get_year(y, areas, v, tf, months, temp_fold)
                for area_name in areas_names:
                        paths[area_name].append(raw_paths[area_name])
                
for i, area_name in enumerate(areas_names):
        print("Merging files for area " + area_name + "...")
        ds=xr.open_mfdataset(paths["NewCaledonia"], parallel=True)
        ds=prep_netcdf(ds)
        delayed = ds.to_netcdf(nc_files[i], compute=False)
        results = delayed.compute(scheduler='threads')
        
shutil.rmtree(temp_fold)
