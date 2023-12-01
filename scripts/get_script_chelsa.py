import fsspec
import xarray as xr

# variables
xmin = 155.86890
xmax = 172.09009
ymin = -22.84806
ymax = -17.39917
variable = "tas"

# test
a = []
month = 1
url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/ncdf/CHELSA_' + \
  variable + '_' + '%02d' % (month,) + '_1981-2010_V.2.1.nc'
fobj = fsspec.open(url)
ds = xr.open_dataset(fobj).chunk({'lat': 500, 'lon': 500})
ds = xr.open_dataset("/home/sschmitt/Téléchargements/CHELSA_tas_01_1981-2010_V.2.1.nc").chunk({'lat': 500, 'lon': 500})
mask_lon = (ds.lon >= xmin) & (ds.lon <= xmax)
mask_lat = (ds.lat >= ymin) & (ds.lat <= ymax)
ds = ds.where(mask_lon & mask_lat, drop = True)
ds.load()
a.append(ds)
    
ds = xr.concat([i for i in a], 'time')

 
 
# get urls
a = []
for month in range(1, 13):
  url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/ncdf/CHELSA_' + \
    variable + '_' + '%02d' % (month,) + '_1981-2010_V.2.1.nc'
    with fsspec.open(url) as fobj:
      ds = xr.open_dataset(fobj).chunk({'lat': 500, 'lon': 500})
      mask_lon = (ds.lon >= xmin) & (ds.lon <= xmax)
      mask_lat = (ds.lat >= ymin) & (ds.lat <= ymax)
      ds = ds.where(mask_lon & mask_lat, drop = True)
      ds.load()
    a.append(ds)
    
ds = xr.concat([i for i in a], 'time')
      
if self.variable_id == "tas" or self.variable_id == 'tasmin' or self.variable_id == 'tasmax':
  res = ds.assign(Band1=ds['Band1'] * 0.1)
if self.variable_id == 'pr':
  res = ds.assign(Band1=ds['Band1'] * 0.1)
  
            



import fsspec
import xarray as xr

class chelsaV2:
    """ 
    Class to download and clip data from the CHELSA V2.1 normals (climatologies)
    for a specific bounding box delimited by minimum and maximum latitude and longitude
    
    :param xmin: Minimum longitude [Decimal degree]
    :param xmax: Maximum longitude [Decimal degree]
    :param ymin: Minimum latitude [Decimal degree]
    :param ymax: Maximum latitude [Decimal degree]
    :param variable_id: id of the variable that needs to be downloaded (e.g. 'tas')

    """
    def __init__(self, xmin, xmax, ymin, ymax, variable_id):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.variable_id = variable_id

    def _crop_ds_(self, ds):
        """
        clip xarray
        
        :param ds: a xarray to_dataset
        :return: clipped xarray
        :rtype: xarray
        """
        mask_lon = (ds.lon >= self.xmin) & (ds.lon <= self.xmax)
        mask_lat = (ds.lat >= self.ymin) & (ds.lat <= self.ymax)
        cropped_ds = ds.where(mask_lon & mask_lat, drop=True)
        return cropped_ds

    def get_chelsa(self):
        """
        download chelsa
        
        :return: cropped xarray
        :rtype: xarray
        """

        a = []
        for month in range(1, 13):
            url = 'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/ncdf/CHELSA_' + self.variable_id + '_' + '%02d' % (
                month,) + '_1981-2010_V.2.1.nc'
            with fsspec.open(url) as fobj:
                ds = xr.open_dataset(fobj).chunk({'lat': 500, 'lon': 500})
                ds = self._crop_ds_(ds)
                ds.load()
            a.append(ds)

        ds = xr.concat([i for i in a], 'time')

        #ds = self._crop_ds_(xr.concat([i for i in a], 'time'))

        # old version using rasterio
        #a = []
        #for month in range(1, 13):
        #    url = 'https://envicloud.os.zhdk.cloud.switch.ch/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/' + self.variable_id + '/CHELSA_' + self.variable_id + '_' + '%02d' % (month,) + '_1981-2010_V.2.1.tif'
        #    a.append(url)

        #ds = self._crop_ds_(xr.concat([xr.open_rasterio(i) for i in a], 'time'))
        if self.variable_id == "tas" or self.variable_id == 'tasmin' or self.variable_id == 'tasmax':
            res = ds.assign(Band1=ds['Band1'] * 0.1)
        if self.variable_id == 'pr':
            res = ds.assign(Band1=ds['Band1'] * 0.1)

        return res
      
chelsa = chelsaV2(155.86890, 172.09009, -22.84806, -17.39917, "tas")

tas = chelsa.get_chelsa()

tas.to_netcdf("test.nc")

tas.to_netcdf(path = "test.nc", mode = "w", format = "NETCDF4", engine = "netcdf4")
