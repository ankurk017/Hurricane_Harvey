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



wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01_no_moving_domain/WRF/test/em_real/'
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d02*"))[::3]


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4_test_may19/'
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d03*"))[::5]


basin = tracks.TrackDataset(basin='north_atlantic',
                            source='ibtracs', include_btk=False)
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

"""
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


"""
# axs.plot(values["Time"], values, "b-")

fig, axs = plt.subplots(2, 1, figsize=(8, 7))

axs[0].plot(ws_timeseries["Time"], ws_timeseries*1.95, "r", label="LULC 2001")

axs[0].plot(obs_ws["Time"],
         obs_ws*1.95, "b", label="LULC 2017")

axs[0].plot(harvey['date'], harvey['vmax'], "k-", label="OBS")
axs[0].plot([harvey['date'][42], harvey['date'][42]], [10, 150], 'b--')
# plt.legend()
axs[0].set_xlabel("Date")
axs[0].set_ylabel("10 m Wind Speed (knots)")
axs[0].set_xlim((harvey['date'][33], harvey['date'][50]))
axs[0].set_ylim((10, 150))


myFmt = mdates.DateFormatter('%m-%d')
myFmt = mdates.DateFormatter('%dZ%H')
axs[0].xaxis.set_major_formatter(myFmt)
#plt.title(ws_error_var)
#plt.tight_layout()
#plt.savefig('../figures/May22_WS.jpeg')

#fig, axs = plt.subplots(figsize=(8, 4))
axs[1].plot(mslp_timeseries["Time"], mslp_timeseries, "r", label="LULC 2001")

axs[1].plot(obs_slp["Time"],
         obs_slp, "b", label="LULC 2017")

axs[1].plot(harvey['date'], harvey['mslp'], "k-", label="OBS")
axs[1].plot([harvey['date'][42], harvey['date'][42]], [880, 1008], 'b--')
# plt.legend()
# axs.set_xticks(rotation=45)
axs[1].set_xlabel("Date")
axs[1].set_ylabel("MSLP (hPa)")
axs[1].set_xlim((harvey['date'][33], harvey['date'][50]))
axs[1].set_ylim((880, 1008))
# myFmt = mdates.DateFormatter('%m-%d')
myFmt = mdates.DateFormatter('%dZ%H')
axs[1].xaxis.set_major_formatter(myFmt)

#plt.title(sp_error_var)
plt.tight_layout()
plt.savefig('../figures/May22_Intensity_nomoving.jpeg')

# plt.show()

# os.system(f'grep -r ATCF {wrf_pruns}/rsl.error.0000  > IDA_track_Moving_Nest.txt')
# wrf_track=pd.read_csv('IDA_track_Moving_Nest.txt',header=None, delimiter=r"\s+")
# plt.plot([datetime.datetime.strptime(dates, '%Y-%m-%d_%H:%M:%S') for dates in wrf_track[1]], wrf_track[4], 'm-')


# PLOT TRACK
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()




fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
coast.plot_coast(ax)
for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')
ax.plot(harvey['lon'], harvey['lat'], 'k-', label='OBS')
ax.plot(track_lon[6:], track_lat[6:], 'r-', label='1 km (9-3-1 km no-moving)')
ax.plot(track_lon_phy, track_lat_phy, 'b-', label='1 km (9-3-1 km moving)')

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([-102, -83])  
ax.set_ylim([18, 35]) 
plt.legend(loc='upper left')

plt.tight_layout()
plt.savefig('../figures/May22_Harvey_track_nomoving.jpeg')
plt.show()

