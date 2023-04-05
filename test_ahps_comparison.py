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

wrfoutfile = sorted(glob.glob(
    '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/wrfout_d01_2017-*26*'))

wrf_pcp = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") +
    getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (
    getvar(Dataset(wrfoutfile[0]), "RAINC") +
    getvar(Dataset(wrfoutfile[0]), "RAINNC")
)

wrf_pcp = sum([getvar(Dataset(wrffile_out), "RAINC") +
              getvar(Dataset(wrffile_out), "RAINNC") for wrffile_out in wrfoutfile[-5:]])


slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)

ahps_fileid = np.where(['20170826' in files for files in ahps_files])[0][0]
ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])


domain_bb = [-100, -91, 24, 33]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500))

fig, ax = plt.subplots(
    1, 2, figsize=(14, 5), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)
cont = ax[0].contourf(ahps_lon, ahps_lat, ahps_pcp,
                      levels=levels, locator=ticker.LogLocator(), cmap='jet')

coast.plot_coast(ax[0])
ax[0].set_xlim((domain_bb[0], domain_bb[1]))
ax[0].set_ylim((domain_bb[2], domain_bb[3]))
ax[0].set_title('AHPS')

cont = ax[1].contourf(wrf_lon, wrf_lat, wrf_pcp, cmap="jet",
                      levels=levels, locator=ticker.LogLocator())

coast.plot_coast(ax[1])
ax[1].set_xlim((domain_bb[0], domain_bb[1]))
ax[1].set_ylim((domain_bb[2], domain_bb[3]))
ax[1].set_title('WRF pre')

cbar = plt.colorbar(cont, ax=ax.ravel(), ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')
plt.savefig('../figures/AHPS_pre_pcp.jpeg')
plt.show()
