from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"

ahps_files = glob.glob(home + "nws_precip*_conus.nc")


ahps_fileid = np.where(['20170827' in files for files in ahps_files])[0][0]
ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])


domain_bb = [-100, -91, 24, 33]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600))

fig, ax = plt.subplots(
    1, 1, figsize=(9, 6), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)
cont = ax.contourf(ahps_lon, ahps_lat, ahps_pcp,
                      levels=levels, locator=ticker.LogLocator(), cmap='jet', extend='both')

coast.plot_coast(ax)
ax.set_xlim((domain_bb[0], domain_bb[1]))
ax.set_ylim((domain_bb[2], domain_bb[3]))
ax.set_title('AHPS')

cbar = plt.colorbar(cont, ax=ax
, ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')
plt.savefig('AHPS_rainfall.jpeg', dpi=400)

#plt.show()




