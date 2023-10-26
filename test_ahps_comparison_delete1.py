from matplotlib import ticker, cm

import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np

import datetime

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

gpm_files = sorted(
    glob.glob("/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/GPM/*nc4")
)[432:650]
gpm_rainfall = xr.open_mfdataset(gpm_files)["precipitationCal"]
gpm_timestep = [
    datetime.datetime(*dates.timetuple()[:6]) for dates in gpm_rainfall["time"].values
]

gpm = xr.open_mfdataset(gpm_files)["precipitationCal"]
gpm_croppped = gpm.isel(time=np.arange(gpm_timestep.index(datetime.datetime(2017, 8, 26, 12, 0)), gpm_timestep.index(datetime.datetime(2017, 8, 27, 12, 0)),1)).sum(dim='time')
#gpm_croppped = gpm.isel(time=np.arange(gpm_timestep.index(datetime.datetime(2017, 8, 27, 0, 0)), gpm_timestep.index(datetime.datetime(2017, 8, 28, 0, 0)),1)).sum(dim='time')

lon = gpm_croppped.lon.values
lat = gpm_croppped.lat.values



domain_bb = [-100, -91, 24, 33]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600))

fig, ax = plt.subplots(
    1, 1, figsize=(9, 6), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)
cont = ax.contourf(lon, lat, gpm_croppped.T,
                      levels=levels, locator=ticker.LogLocator(), cmap='jet', extend='both')

coast.plot_coast(ax)
ax.set_xlim((domain_bb[0], domain_bb[1]))
ax.set_ylim((domain_bb[2], domain_bb[3]))
ax.set_title('GPM')

cbar = plt.colorbar(cont, ax=ax
, ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')

plt.savefig('GPM_2612-2712.jpeg', dpi=400)
#plt.show()




