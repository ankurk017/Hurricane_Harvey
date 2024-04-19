from src.wrf_src import wrf_assign_coords
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


wrfoutfile = sorted(
    glob.glob(
        "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/wrfout_d01_2017-*"
    )
)
wrf_pcp = getvar(Dataset(wrfoutfile[-1]), "RAINC") + getvar(
    Dataset(wrfoutfile[-1]), "RAINNC"
)


wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location

var_timeseries_pre = []
var_timeseries_post = []

for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_pre))):

    wrf_ncfile_pre = Dataset(wrfoutfile_pre[wrf_files_id])
    var_pre = wrf_assign_coords(
        getvar(wrf_ncfile_pre, "RAINC") + getvar(wrf_ncfile_pre, "RAINNC")
    )

    lat_id = np.where(np.logical_and(
        var_pre['south_north'].values > location[1]-0.1, var_pre['south_north'].values <= location[1]+0.1))[0]
    lon_id = np.where(np.logical_and(
        var_pre['west_east'].values > location[0]-0.1, var_pre['west_east'].values <= location[0]+0.1))[0]

    var_timeseries_pre.append(var_pre.isel(
        south_north=lat_id, west_east=lon_id).mean())

    wrf_ncfile_post = Dataset(wrfoutfile_post[wrf_files_id])
    var_post = wrf_assign_coords(
        getvar(wrf_ncfile_post, "RAINC") + getvar(wrf_ncfile_post, "RAINNC")
    )
    var_timeseries_post.append(var_post.isel(
        south_north=lat_id, west_east=lon_id).mean())

var_timeseries_pre_merged = xr.concat(var_timeseries_pre, dim="time")
var_timeseries_post_merged = xr.concat(var_timeseries_post, dim="time")

####### GPM ####
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
# axs.plot([gpm_timestep[102], gpm_timestep[102]], [-10, 520], 'k--')
axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm/hr)")
plt.legend()
plt.xticks(rotation=25)
plt.xlim((var_timeseries_post_merged["Time"]
         [0], var_timeseries_post_merged["Time"][-1]))
plt.ylim((0, 510))

plt.tight_layout()
# plt.savefig('../figures/WRF_GPM_rainfall.jpeg')
plt.show()
