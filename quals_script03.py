from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from wrf import getvar
import glob
import xarray as xr
import progressbar
from src.wrf_src import wrf_assign_coords
from src.coast import plot_coast
import numpy as np

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[48]

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[48]

var_name = "slp"
location = (-97.061, 27.8339)  # Landfall location


wrf_ncfile = Dataset(wrfoutfile_pre)
var = wrf_assign_coords(getvar(wrf_ncfile, 'uvmet10'))
u = var.sel(u_v='u')
v = var.sel(u_v='v')
ws_pre = np.sqrt(u**2 + v**2)


wrf_ncfile = Dataset(wrfoutfile_post)
var = wrf_assign_coords(getvar(wrf_ncfile, 'uvmet10'))
u = var.sel(u_v='u')
v = var.sel(u_v='v')
ws_post = np.sqrt(u**2 + v**2)

diff = ws_post-ws_pre


diff.plot()
plt.scatter(-96.9, 28, color='b')
# plt.xlim((-105, -85))

# plt.ylim((22, 30))
plt.show()


# ax = plt.axes(projection=ccrs.PlateCarree())
# diff.plot.contourf(ax=ax, transform=ccrs.PlateCarree())
# plot_coast(ax)
# plt.show()
