from self_utils import coast
from netCDF4 import Dataset
import matplotlib.dates as mdates
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import wrf
import pandas as pd
from wrf import (
    getvar,
    interplevel,
    to_np,
    latlon_coords,
    get_cartopy,
    cartopy_xlim,
    cartopy_ylim,
    ALL_TIMES,
)
import glob
import pandas as pd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import re
import os
import xarray as xr
import datetime
import progressbar

import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

import tropycal.tracks as tracks


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/'
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[25:25+48]

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4_test_may19/'
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

# wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy02/'

# wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01/'
# wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[25:25+48]




import numpy as np

def calculate_pressure_location(sea_level_pressure):
    
    min_pressure_coords = np.unravel_index(np.argmin(sea_level_pressure), sea_level_pressure.shape)

    pressure_sum = np.sum(sea_level_pressure)
    centroid_coords = np.array([0, 0])
    for i in range(sea_level_pressure.shape[0]):
        for j in range(sea_level_pressure.shape[1]):
            centroid_coords[0] += i * sea_level_pressure[i][j]
            centroid_coords[1] += j * sea_level_pressure[i][j]
    centroid_coords1 = centroid_coords / pressure_sum

    #return np.round(centroid_coords1).astype(int) 
    return centroid_coords1 



basin = tracks.TrackDataset(basin='north_atlantic',
                            source='ibtracs', include_btk=False)
harvey = basin.get_storm(('harvey', 2017))


var_name = "slp"


track_lon = []
track_lat = []
track_lon_centroid = []
track_lat_centroid = []
slp_min = []
ws_max = []
for timeid in progressbar.progressbar(range(len(wrfoutfile_pre))):
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile_pre[timeid])
    slp = getvar(wrf_ncfile, var_name)
    slp_min.append(slp.min())
    ws = getvar(wrf_ncfile, 'uvmet10_wspd_wdir').sel(wspd_wdir='wspd')
    ws_max.append(ws.max())
    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))

    track_lon.append(lons[np.where(slp == slp.min())].values[0][0])
    track_lat.append(lats[np.where(slp == slp.min())].values[0][0])

#    track_lon_centroid.append(lons[calculate_pressure_location(slp.values)].values[0][0])
#    track_lat_centroid.append(lats[calculate_pressure_location(slp.values)].values[0][0])

mslp_timeseries = xr.concat(slp_min, dim="time")
ws_timeseries = xr.concat(ws_max, dim="time")


# PLOT TRACK

fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
coast.plot_coast(ax)
ax.plot(harvey['lon'], harvey['lat'], 'k-', label='OBS')
ax.plot(track_lon[6:], track_lat[6:], 'r-', label='2001')
#ax.plot(track_lon_centroid[6:], track_lat_centroid[6:], 'b-o', label='centroid')
#ax.plot(track_lon_phy[6:-3], track_lat_phy[6:-3], 'b-', label='2017')
# plt.legend()
plt.tight_layout()
plt.show()











