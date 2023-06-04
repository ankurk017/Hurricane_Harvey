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


wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location

box = 0.2
var_timeseries_pre = []
var_timeseries_post = []

for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_pre) - 1)):

    wrf_ncfile_pre = Dataset(wrfoutfile_pre[wrf_files_id + 1])
    wrf_ncfile_pre_tminus1 = Dataset(wrfoutfile_pre[wrf_files_id])

    var_pre_tmp = wrf_assign_coords(
        getvar(wrf_ncfile_pre, "RAINC") + getvar(wrf_ncfile_pre, "RAINNC")
    ) - wrf_assign_coords(
        getvar(wrf_ncfile_pre_tminus1, "RAINC")
        + getvar(wrf_ncfile_pre_tminus1, "RAINNC")
    )

    var_pre = wrf_assign_coords(getvar(wrf_ncfile_pre, "RAINC"))
    var_pre.values = var_pre_tmp.values

    lat_id = np.where(
        np.logical_and(
            var_pre["south_north"].values > location[1] - box,
            var_pre["south_north"].values <= location[1] + box,
        )
    )[0]
    lon_id = np.where(
        np.logical_and(
            var_pre["west_east"].values > location[0] - box,
            var_pre["west_east"].values <= location[0] + box,
        )
    )[0]

    var_timeseries_pre.append(var_pre.isel(south_north=lat_id, west_east=lon_id).mean())

    wrf_ncfile_post = Dataset(wrfoutfile_post[wrf_files_id + 1])
    wrf_ncfile_post_tminus1 = Dataset(wrfoutfile_post[wrf_files_id])
    var_post_tmp = wrf_assign_coords(
        getvar(wrf_ncfile_post, "RAINC") + getvar(wrf_ncfile_post, "RAINNC")
    ) - wrf_assign_coords(
        getvar(wrf_ncfile_post_tminus1, "RAINC")
        + getvar(wrf_ncfile_post_tminus1, "RAINNC")
    )

    var_post = wrf_assign_coords(getvar(wrf_ncfile_post, "RAINC"))
    var_post.values = var_post_tmp.values

    var_timeseries_post.append(
        var_post.isel(south_north=lat_id, west_east=lon_id).mean()
    )

var_timeseries_pre_merged = xr.concat(var_timeseries_pre, dim="Time")
var_timeseries_post_merged = xr.concat(var_timeseries_post, dim="Time")

####### GPM ####
gpm_files = sorted(
    glob.glob("/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/GPM/*nc4")
)[432:650]
gpm_rainfall = xr.open_mfdataset(gpm_files)["precipitationCal"].interp(
    lon=location[0], lat=location[1]
)
gpm_timestep = [
    datetime.datetime(*dates.timetuple()[:6]) for dates in gpm_rainfall["time"].values
]
gpm = xr.open_mfdataset(gpm_files)["precipitationCal"]
gpm_rainfall_region = (
    gpm.isel(
        lon=np.where(
            np.logical_and(
                gpm["lon"] > location[0] - box, gpm["lon"] <= location[0] + box
            )
        )[0],
        lat=np.where(
            np.logical_and(
                gpm["lat"] > location[1] - box, gpm["lat"] <= location[1] + box
            )
        )[0],
    )
    .mean(dim="lon")
    .mean(dim="lat")
)

#################

fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(gpm_timestep, np.cumsum(gpm_rainfall_region.values), "k-", label="GPM")
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
# axs.plot([gpm_timestep[102], gpm_timestep[102]], [-10, 520], 'k--')
axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm/hr)")
plt.legend()
plt.xticks(rotation=25)
plt.xlim(
    (var_timeseries_post_merged["Time"][0], var_timeseries_post_merged["Time"][-1])
)
# plt.ylim((0, 510))
plt.tight_layout()
# plt.savefig('../figures/WRF_GPM_rainfall.jpeg')


fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(gpm_timestep, gpm_rainfall_region.values, "k-", label="GPM")
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
plt.xlim(
    (var_timeseries_post_merged["Time"][0], var_timeseries_post_merged["Time"][-1])
)
# plt.ylim((0, 510))
plt.tight_layout()
# plt.savefig('../figures/WRF_GPM_rainfall.jpeg')
plt.show()





