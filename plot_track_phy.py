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

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_GFS/LULC_2001/WRF/test/em_real/'
wrf_pruns = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_NAM/LULC_2001/WRF/test/em_real/'

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_pre/WRF/WRFV4/'
wrf_pruns = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_post/WRF/WRFV4/'


wrfoutfile = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrfoutfile_physics = sorted(glob.glob(wrf_pruns + "wrfout_d02*"))



basin = tracks.TrackDataset(basin='north_atlantic',source='hurdat',include_btk=False)

harvey = basin.get_storm(('harvey',2017))



wrf_ncfile = Dataset(wrfoutfile[0])

time_id = 30
var_name = "slp"


track_lon = []
track_lat = []
slp_min = []
ws_max = []
for timeid in progressbar.progressbar(range(len(wrfoutfile))):
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile[timeid])
    slp = getvar(wrf_ncfile, var_name)
    slp_min.append(slp.min())
    ws = getvar(wrf_ncfile, 'uvmet10_wspd_wdir').sel(wspd_wdir='wspd')
    ws_max.append(ws.max())
    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
    track_lon.append(lons[np.where(slp == slp.min())].values[0][0])
    track_lat.append(lats[np.where(slp == slp.min())].values[0][0])

mslp_timeseries = xr.concat(slp_min, dim="time")
ws_timeseries = xr.concat(ws_max, dim="time")
#wrf_ida_lonlat_land = [is_land(lon, lat) for lon, lat in zip(track_lon, track_lat)]
#ida_wrf_landfall_time = mslp_timeseries.isel(time=42)["Time"].values


# PHYSICS
ptrack_lon = []
ptrack_lat = []
slp_min = []
ws_max = []
for timeid in progressbar.progressbar(range(len(wrfoutfile_physics))):
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile_physics[timeid])
    slp = getvar(wrf_ncfile, var_name)
    slp_min.append(slp.min())
    ws = getvar(wrf_ncfile, 'uvmet10_wspd_wdir').sel(wspd_wdir='wspd')
    ws_max.append(ws.max())
    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
    ptrack_lon.append(lons[np.where(slp == slp.min())].values[0][0])
    ptrack_lat.append(lats[np.where(slp == slp.min())].values[0][0])

mslp_timeseries_physics = xr.concat(slp_min, dim="time")
ws_timeseries_physics = xr.concat(ws_max, dim="time")


# PHYSICS
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)

plt.plot(track_lon, track_lat, "r", label="2001")
plt.plot(ptrack_lon, ptrack_lat, "b", label="2020")

plt.plot(harvey['lon'], harvey['lat'], "k",
         transform=ccrs.PlateCarree(), label="OBS")

plt.legend()

plt.grid(True)
ax.set_ylim(10, 50)
ax.set_xlim(-100, -60)
gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=True,
    linewidth=2,
    color="gray",
    alpha=0.5,
    linestyle="--",
)
gl.top_labels = True
gl.left_ylabels = True
gl.xlines = True
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
# plt.savefig('../figures/Track_comparison_corral.jpeg')


fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(ws_timeseries["Time"], ws_timeseries*1.95, "r", label="2001")
axs.plot(ws_timeseries_physics["Time"],
         ws_timeseries_physics*1.95, "b", label="2020")
axs.plot(harvey['date'], harvey['vmax'], "k-", label="OBS")
plt.legend()
axs.set_xlabel("Date")
axs.set_ylabel("10 m Wind Speed (knots)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
myFmt = mdates.DateFormatter('%m-%d')
axs.xaxis.set_major_formatter(myFmt)
plt.tight_layout()

#plt.savefig('../figures/Oct05_Track_WS_comparison_corral.jpeg')

fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(mslp_timeseries["Time"], mslp_timeseries, "r", label="2001")
axs.plot(mslp_timeseries_physics["Time"],
         mslp_timeseries_physics, "b", label="2020")
axs.plot(harvey['date'], harvey['mslp'], "k-", label="OBS")
plt.legend()
#axs.set_xticks(rotation=45)
axs.set_xlabel("Date")
axs.set_ylabel("MSLP (hPa)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
plt.tight_layout()
myFmt = mdates.DateFormatter('%m-%d')
axs.xaxis.set_major_formatter(myFmt)

#plt.savefig('../figures/Oct05_Track_MSLP_comparison_corral.jpeg')
plt.show()

#os.system(f'grep -r ATCF {wrf_pruns}/rsl.error.0000  > IDA_track_Moving_Nest.txt')
#wrf_track=pd.read_csv('IDA_track_Moving_Nest.txt',header=None, delimiter=r"\s+")
#plt.plot([datetime.datetime.strptime(dates, '%Y-%m-%d_%H:%M:%S') for dates in wrf_track[1]], wrf_track[4], 'm-')



plt.show()
