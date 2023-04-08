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

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/'
wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_Quals/WRF_cropland/WRF/'
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d02*"))[::3]


basin = tracks.TrackDataset(basin='north_atlantic',
                            source='hurdat', include_btk=False)
harvey = basin.get_storm(('harvey', 2017))


var_name = "slp"


track_lon = []
track_lat = []
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

mslp_timeseries = xr.concat(slp_min, dim="time")
ws_timeseries = xr.concat(ws_max, dim="time")
# wrf_ida_lonlat_land = [is_land(lon, lat) for lon, lat in zip(track_lon, track_lat)]
# ida_wrf_landfall_time = mslp_timeseries.isel(time=42)["Time"].values


track_lon_phy = []
track_lat_phy = []
slp_min = []
ws_max = []
for timeid in progressbar.progressbar(range(len(wrfoutfile_post))):
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile_post[timeid])
    slp = getvar(wrf_ncfile, var_name)
    slp_min.append(slp.min())
    ws = getvar(wrf_ncfile, 'uvmet10_wspd_wdir').sel(wspd_wdir='wspd')
    ws_max.append(ws.max())
    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
    track_lon_phy.append(lons[np.where(slp == slp.min())].values[0][0])
    track_lat_phy.append(lats[np.where(slp == slp.min())].values[0][0])

obs_slp = xr.concat(slp_min, dim="time")
obs_ws = xr.concat(ws_max, dim="time")
# wrf_ida_lonlat_land = [is_land(lon, lat) for lon, lat in zip(track_lon, track_lat)]
# ida_wrf_landfall_time = mslp_timeseries.isel(time=42)["Time"].values

# calculating the time series error

model_ws = np.array([harvey['vmax'][np.where(harvey['date'] == pd.to_datetime(
    str(val.values)).to_pydatetime())[0]] for val in obs_slp["Time"]], dtype='object')
model2_ws_error = (model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]-(
    obs_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]*1.95).values).mean()

model_ws = np.array([harvey['vmax'][np.where(harvey['date'] == pd.to_datetime(str(
    val.values)).to_pydatetime())[0]] for val in mslp_timeseries["Time"]], dtype='object')
model1_ws_error = (model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]-(
    ws_timeseries[np.array([ws.shape[0] != 0 for ws in model_ws])]*1.95).values).mean()

ws_error_var = f'WS Error: {model1_ws_error}, {model2_ws_error}'


model_ws = np.array([harvey['mslp'][np.where(harvey['date'] == pd.to_datetime(
    str(val.values)).to_pydatetime())[0]] for val in obs_slp["Time"]], dtype='object')
model2_sp_error = (model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]-(
    obs_slp[np.array([ws.shape[0] != 0 for ws in model_ws])]).values).mean()

model_ws = np.array([harvey['mslp'][np.where(harvey['date'] == pd.to_datetime(str(
    val.values)).to_pydatetime())[0]] for val in mslp_timeseries["Time"]], dtype='object')
model1_sp_error = (model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]-(
    mslp_timeseries[np.array([ws.shape[0] != 0 for ws in model_ws])]).values).mean()

sp_error_var = f'MSLP Error: {model1_sp_error}, {model2_sp_error}'


# axs.plot(values["Time"], values, "b-")

fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(ws_timeseries["Time"], ws_timeseries*1.95, "r", label="LULC 2001")

axs.plot(obs_ws["Time"],
         obs_ws*1.95, "b", label="LULC 2017")

axs.plot(harvey['date'], harvey['vmax'], "k-", label="OBS")
# plt.legend()
axs.set_xlabel("Date")
axs.set_ylabel("10 m Wind Speed (knots)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
myFmt = mdates.DateFormatter('%m-%d')
myFmt = mdates.DateFormatter('%dZ%H')
axs.xaxis.set_major_formatter(myFmt)
plt.title(ws_error_var)
plt.tight_layout()
# plt.savefig('../figures/WS.jpeg')

fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(mslp_timeseries["Time"], mslp_timeseries, "r", label="LULC 2001")

axs.plot(obs_slp["Time"],
         obs_slp, "b", label="LULC 2017")

axs.plot(harvey['date'], harvey['mslp'], "k-", label="OBS")
# plt.legend()
# axs.set_xticks(rotation=45)
axs.set_xlabel("Date")
axs.set_ylabel("MSLP (hPa)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
# myFmt = mdates.DateFormatter('%m-%d')
myFmt = mdates.DateFormatter('%dZ%H')
axs.xaxis.set_major_formatter(myFmt)

plt.title(sp_error_var)
plt.tight_layout()
# plt.savefig('../figures/MSLP.jpeg')

# plt.show()

# os.system(f'grep -r ATCF {wrf_pruns}/rsl.error.0000  > IDA_track_Moving_Nest.txt')
# wrf_track=pd.read_csv('IDA_track_Moving_Nest.txt',header=None, delimiter=r"\s+")
# plt.plot([datetime.datetime.strptime(dates, '%Y-%m-%d_%H:%M:%S') for dates in wrf_track[1]], wrf_track[4], 'm-')


# PLOT TRACK


fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
coast.plot_coast(ax)
ax.plot(harvey['lon'], harvey['lat'], 'k-', label='OBS')
ax.plot(track_lon[6:], track_lat[6:], 'r-', label='2001')
ax.plot(track_lon_phy[6:-3], track_lat_phy[6:-3], 'b-', label='2017')
# plt.legend()
plt.tight_layout()
plt.show()
