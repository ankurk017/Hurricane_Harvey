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

from netCDF4 import Dataset
import matplotlib.pyplot as plt

from wrf import getvar
import glob
import xarray as xr
import progressbar
from src.wrf_src import wrf_assign_coords


wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location

var_timeseries_pre = []
var_timeseries_post = []

for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile))):

    wrf_ncfile_pre = Dataset(wrfoutfile_pre[wrf_files_id])
    var_pre = wrf_assign_coords(
        getvar(wrf_ncfile_pre, "RAINC") + getvar(wrf_ncfile_pre, "RAINNC")
    )
    var_timeseries_pre.append(
        var_pre.interp(south_north=location[1], west_east=location[0], method="nearest")
    )

    wrf_ncfile_post = Dataset(wrfoutfile_post[wrf_files_id])
    var_post = wrf_assign_coords(
        getvar(wrf_ncfile_post, "RAINC") + getvar(wrf_ncfile_post, "RAINNC")
    )
    var_timeseries_post.append(
        var_post.interp(
            south_north=location[1], west_east=location[0], method="nearest"
        )
    )

var_timeseries_pre_merged = xr.concat(var_timeseries_pre, dim="time")
var_timeseries_post_merged = xr.concat(var_timeseries_post, dim="time")

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
#################

fig, axs = plt.subplots(1, 1, figsize=(7, 4))
axs.plot(gpm_timestep, np.cumsum(gpm_rainfall.values), 'k-', label='GPM')
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
plt.tight_layout()
plt.show()
