from src.wrf_src import wrf_assign_coords
import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords, interplevel
from self_utils import ahps, coast
import numpy as np

import datetime

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})



wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = wrf_runs + "wrfout_d02_2017-08-24_18:00:00"

var_name = "slp"


ncfile = Dataset(wrfoutfile_pre)
var_pre = wrf_assign_coords( getvar(ncfile, "wa")  )


# Extract the pressure, geopotential height, and wind variables
p = getvar(ncfile, "pressure")
z = getvar(ncfile, "z", units="dm")

wa = getvar(ncfile, "wa", units="kt")
wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]


fig, axs = plt.subplots(1, 2, figsize=(15, 6))
interplevel(wspd, p, 200).plot(ax=axs[0])
interplevel(wa, p, 200).plot(ax=axs[1])

fig, axs = plt.subplots(1, 2, figsize=(15, 6))
interplevel(wspd, p, 850).plot(ax=axs[0])
interplevel(wa, p, 850).plot(ax=axs[1])
plt.show()


