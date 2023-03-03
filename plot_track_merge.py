from netCDF4 import Dataset
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



plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


wrfoutfile = sorted(glob.glob("/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/LULC_2001/WRF/test/em_real/wrfout_d02*"))

wrf_ncfile = Dataset(wrfoutfile[0])

time_id = 30
var_name = "slp"


track_lon1 = []
track_lat1 = []
slp_min = []
ws_max = []
for timeid in progressbar.progressbar(range(len(wrfoutfile))):
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile[timeid])
    slp = getvar(wrf_ncfile, "slp")
    slp_min.append(slp.min())
    ws = getvar(wrf_ncfile, 'uvmet10_wspd_wdir').sel(wspd_wdir='wspd')
    ws_max.append(ws.max())
    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
    track_lon1.append(lons[np.where(slp == slp.min())].values[0][0])
    track_lat1.append(lats[np.where(slp == slp.min())].values[0][0])

mslp_timeseries = xr.concat(slp_min, dim="time")
ws_timeseries = xr.concat(ws_max, dim="time")


# PHYSICS
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)

plt.plot(track_lon1, track_lat1, "b", label="Sim 1")

#plt.plot(-np.array(ida_lon), ida_lat, "k",
#         transform=ccrs.PlateCarree(), label="OBS")
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
#plt.savefig('../figures/L0_P1.jpeg')
plt.show()
