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
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d02*"))[25:25+48]

#wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy02/'

#wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01/'
#wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[25:25+48]


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01/'
wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/'
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[25:25+48]

#wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01c_new_config/'
#wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01d_new_config_from_paper/'
#wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/test_phy/WRF_9-3-51/WRFV4_phy01_no_moving_domain/WRF/test/em_real/'
#wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d02*"))



basin = tracks.TrackDataset(basin='north_atlantic',source='hurdat',include_btk=False)
harvey = basin.get_storm(('harvey',2017))

var_name = "slp"

slp_min = []
tilt = []
for timeid in progressbar.progressbar(range(len(wrfoutfile_pre))):
    
    #    print(timeid)
    wrf_ncfile = Dataset(wrfoutfile_pre[timeid])
    slp = getvar(wrf_ncfile, var_name)
    slp = getvar(wrf_ncfile, var_name)
    slp_min.append(slp.min())

    p = getvar(wrf_ncfile, "pressure")
    ua = getvar(wrf_ncfile, "ua", units="kt")
    va = getvar(wrf_ncfile, "va", units="kt")
    ws = np.sqrt(ua**2 + va**2)
    lats, lons = latlon_coords(va)

    cen_1000 = (lons[np.where(interplevel(ws, p, 1000)==interplevel(ws, p, 1000).min())].values, lats[np.where(interplevel(ws, p, 1000)==interplevel(ws, p, 1000).min())].values) 
    cen_200 = (lons[np.where(interplevel(ws, p, 200)==interplevel(ws, p, 200).min())].values, lats[np.where(interplevel(ws, p, 200)==interplevel(ws, p, 200).min())].values) 

    tilt.append(np.sqrt((cen_1000[0]-cen_200[0])**2  + (cen_1000[1]-cen_200[1])**2))

obs_slp = xr.concat(slp_min, dim="time")

fig, axs = plt.subplots(figsize=(10, 5))

axs.plot(harvey['date'], harvey['mslp'], "k-", label="OBS")
axs.set_xlabel("Date (DDZHH)")
axs.set_ylabel("MSLP (hPa)")
axs.set_xlim((harvey['date'][33], harvey['date'][50]))

#myFmt = mdates.DateFormatter('%m-%d')
myFmt = mdates.DateFormatter('%dZ%H')
axs.xaxis.set_major_formatter(myFmt)
#axs.set_yticks(np.arange(0, 150, 30))

color = 'tab:red'
ax2  = axs.twinx()
#ax2.plot(obs_slp['Time'], np.array(tilt).squeeze())
ax2.scatter(obs_slp['Time'], pd.DataFrame(np.array(tilt).squeeze()).rolling(window=3, min_periods=1).mean().values, color=color, marker='o')
axs.grid(color = 'black', linestyle = '--', linewidth = 0.5)
ax2.set_ylabel(r'Vortex tilt ($^o$)', color=color)
ax2.tick_params(axis='y', labelcolor=color)

ax2.xaxis.set_major_formatter(myFmt)
plt.tight_layout()
plt.savefig('../figures/Hurricane_Vortex_tilt.jpeg')
plt.show()

