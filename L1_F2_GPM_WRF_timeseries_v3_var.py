from src.wrf_src import wrf_assign_coords, plot_bb, area_mean
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



home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre_UCM/WRF/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post_UCM/WRF/test/em_real/wrfout_d02_2017-*'))



var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location
urban_change = {"south_north": 29.74, "west_east": -95.49, "box":1.2}

plot_bb(wrfoutfile_pre[12], location=urban_change)
plt.savefig('../figures/rainfall/domain.jpeg')


var_timeseries_pre = []
var_timeseries_post = []

for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_pre))):

    wrf_ncfile_pre = Dataset(wrfoutfile_pre[wrf_files_id])
    var_pre = wrf_assign_coords(
        getvar(wrf_ncfile_pre, "wspd_wdir10").sel(wspd_wdir='wspd') 
    )
    wrf_ncfile_post = Dataset(wrfoutfile_post[wrf_files_id])
    var_post = wrf_assign_coords(
        getvar(wrf_ncfile_post, "wspd_wdir10").sel(wspd_wdir='wspd') 
    )


    var_timeseries_pre.append(area_mean(var_pre, location=urban_change))
    var_timeseries_post.append(area_mean(var_post, location=urban_change))

var_timeseries_pre_merged = xr.concat(var_timeseries_pre, dim="time")
var_timeseries_post_merged = xr.concat(var_timeseries_post, dim="time")


#################

fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(
    var_timeseries_pre_merged["Time"],
    var_timeseries_pre_merged,
    "b-",
    label="LULC 2001",
)
axs.plot(
    var_timeseries_post_merged["Time"],
    var_timeseries_post_merged,
    "r-",
    label="LULC 2017",
)
axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm/hr)")
plt.legend()
plt.xticks(rotation=25)
#plt.xlim((var_timeseries_post_merged["Time"]
#         [0], var_timeseries_post_merged["Time"][-1]))
plt.tight_layout()


fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(
    var_timeseries_pre_merged["Time"],
    np.cumsum(var_timeseries_pre_merged),
    "b-",
    label="LULC 2001",
)
axs.plot(
    var_timeseries_post_merged["Time"],
    np.cumsum(var_timeseries_post_merged),
    "r-",
    label="LULC 2017",
)
axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm/hr)")
plt.legend()
plt.xticks(rotation=25)
plt.xlim((var_timeseries_post_merged["Time"]
         [0], var_timeseries_post_merged["Time"][-1]))
plt.tight_layout()

plt.show()





