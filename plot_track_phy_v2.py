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


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_pre/WRF/WRFV4/'
wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_pre/WRF_9-3-1-51_levels/WRFV4/wrfout/'
wrf_runs = '/nas/rstor/akumar/India/SST-2_runs/2020_AMPHAN_test_for_Harvey/WRF/WRFV4_3doms/wrf/'

# wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_pre/WRF_new_not_moving/WRFV4/'


wrfoutfile = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

basin = tracks.TrackDataset(basin='north_atlantic',
                            source='hurdat', include_btk=False)

harvey = basin.get_storm(('harvey', 2017))


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
# wrf_ida_lonlat_land = [is_land(lon, lat) for lon, lat in zip(track_lon, track_lat)]
# ida_wrf_landfall_time = mslp_timeseries.isel(time=42)["Time"].values

wrffile = "/nas/rstor/akumar/India/SST-2_runs/2020_AMPHAN_test_for_Harvey/WRF/WRFV4/wrf/wrfout_d02_2017-08-24_00:00:00"

ncfile = Dataset(wrffile)

total_time_ids = wrf.extract_times(ncfile, timeidx=ALL_TIMES).shape[0]
obs_slp = xr.concat(
    [getvar(ncfile, "slp", timeidx=timeid).min()
     for timeid in np.arange(0, total_time_ids, 1)],
    dim="Time",
)
obs_ws = xr.concat(
    [getvar(ncfile, "uvmet10_wspd_wdir", timeidx=timeid).sel(
        wspd_wdir='wspd').max() for timeid in np.arange(0, total_time_ids, 1)],
    dim="Time",
)

# axs.plot(values["Time"], values, "b-")

fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(ws_timeseries["Time"], ws_timeseries*1.95, "r", label="WRF_CNTL")

axs.plot(obs_ws["Time"],
         obs_ws*1.95, "b", label="WRF_PHY")

axs.plot(harvey['date'], harvey['vmax'], "k-", label="OBS")
plt.legend()
axs.set_xlabel("Date")
axs.set_ylabel("10 m Wind Speed (knots)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
myFmt = mdates.DateFormatter('%m-%d')
axs.xaxis.set_major_formatter(myFmt)
plt.tight_layout()

# plt.savefig('../figures/Oct05_Track_WS_comparison_corral.jpeg')

fig, axs = plt.subplots(figsize=(10, 5))
axs.plot(mslp_timeseries["Time"], mslp_timeseries, "r", label="WRF_CNTL")

axs.plot(obs_slp["Time"],
         obs_slp, "b", label="WRF_PHY")

axs.plot(harvey['date'], harvey['mslp'], "k-", label="OBS")
plt.legend()
# axs.set_xticks(rotation=45)
axs.set_xlabel("Date")
axs.set_ylabel("MSLP (hPa)")
axs.set_xlim((harvey['date'][33], harvey['date'][54]))
plt.tight_layout()
myFmt = mdates.DateFormatter('%m-%d')
axs.xaxis.set_major_formatter(myFmt)

# plt.savefig('../figures/Oct05_Track_MSLP_comparison_corral.jpeg')
plt.show()

# os.system(f'grep -r ATCF {wrf_pruns}/rsl.error.0000  > IDA_track_Moving_Nest.txt')
# wrf_track=pd.read_csv('IDA_track_Moving_Nest.txt',header=None, delimiter=r"\s+")
# plt.plot([datetime.datetime.strptime(dates, '%Y-%m-%d_%H:%M:%S') for dates in wrf_track[1]], wrf_track[4], 'm-')
