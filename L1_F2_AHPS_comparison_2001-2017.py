import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"


wrfoutfile = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/wrfout_d01_2017-*26*'))
wrf_pcp_pre = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") + getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (
    getvar(Dataset(wrfoutfile[0]), "RAINC") + getvar(Dataset(wrfoutfile[0]), "RAINNC")
)

wrfoutfile = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/wrfout_d01_2017-*26*'))
wrf_pcp_post = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") + getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (
    getvar(Dataset(wrfoutfile[0]), "RAINC") + getvar(Dataset(wrfoutfile[0]), "RAINNC")
)




slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)


from matplotlib import ticker, cm


domain_bb = [-100, -91, 24, 33]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500))

fig, ax = plt.subplots(
    1, 2, figsize=(14, 5), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)

cont = ax[1].contourf(wrf_lon, wrf_lat, wrf_pcp_post-wrf_pcp_pre,  cmap="jet", levels=[-200, -150, -100, -50, 0, 50, 100, 150, 200])

coast.plot_coast(ax[1])
ax[1].set_xlim((domain_bb[0], domain_bb[1]))
ax[1].set_ylim((domain_bb[2], domain_bb[3]))
ax[1].set_title('WRF (LULC 2017)')

cbar = plt.colorbar(cont, ax=ax.ravel(), ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')
plt.savefig('../figures/AHPS_pre-post_pcp.jpeg')
plt.show()
"""






