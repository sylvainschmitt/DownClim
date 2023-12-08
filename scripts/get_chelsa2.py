# log
log_file = open(snakemake.log[0],"w")
sys.stderr = sys.stdout = log_file

# variables
area_file = snakemake.input[0]
area = snakemake.params.area
variables = snakemake.params.variables
time_frequency = snakemake.params.time_frequency
base_years = snakemake.params.base_years
nc_file = snakemake.output[0]

# test
# area_file = "results/countries/New-Caledonia.shp"
# area="New-Caledonia"
# variables = ["tas", "pr"]
# time_frequency = "mon"
# base_years = "1980-1981"
# nc_file = "test.nc"

# libs
import xarray as xr
import rioxarray as rio
import pandas as pd
import geopandas
import datetime

# function
def get_year(year, area, var, time_freq, month_max):
       mm = month_max + 1
       a = []
       for month in range(1, mm):
                url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/' + time_freq + '/' + var + "/" \
                        + 'CHELSA_' + var + '_' + '%02d' % (month,) + '_' + str(year) + "_V.2.1.tif"
                ds = rio.open_rasterio(url, decode_coords="all").to_dataset('band').rename_vars({1 : "var"})
                mask_lon = (ds.x >= area.bounds.minx[0]) & (ds.x <= area.bounds.maxx[0])
                mask_lat = (ds.y >= area.bounds.miny[0]) & (ds.y <= area.bounds.maxy[0])
                ds = ds.where(mask_lon & mask_lat, drop = True)
                a.append(ds)
       ds_year = xr.concat([i for i in a], pd.Index(pd.date_range(datetime.datetime(year, 1, 1), periods=month_max, freq="M"), name="time"))
       ds_year = ds_year[["time", "x", "y", "var"]]
       return(ds_year)

def get_var(years, area, var, time_freq, month_max):
        a = []
        for year in years:
                ds = get_year(year, area, var, time_freq, month_max)
                a.append(ds)
        ds_var = xr.concat(a, "time")
        ds_var["var"] = ds_var["var"] * 0.1
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
        return(ds_var)

# code
area = geopandas.read_file(area_file)

base_years = list(map(int, base_years.split("-")))
years = list(range(base_years[0], base_years[1]+1))

if time_frequency == "mon":
        tf = "monthly"
        months = 12
        # months = 2 # for tests
if time_frequency != "mon":
        raise Exception("Currently only monthly time frequency available!")

a = []
for var in variables:
        ds = get_var(years, area, var, tf, months)
        a.append(ds)
ds_all = xr.merge(a)

# from matplotlib import pyplot as plt # check with plot
# ds_all.pr[0].plot()
# plt.show()

ds_all.to_netcdf(nc_file)
