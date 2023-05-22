from netCDF4 import Dataset
import matplotlib.pyplot as plt

from wrf import getvar
import glob
import xarray as xr
import progressbar
from src.wrf_src import wrf_assign_coords

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location


var_timeseries = []
for wrf_files in progressbar.progressbar(wrfoutfile_pre):
    wrf_ncfile = Dataset(wrf_files)
    var = wrf_assign_coords(getvar(wrf_ncfile, var_name))
    var_timeseries.append(
        var.interp(south_north=location[1],
                   west_east=location[0], method="nearest")
    )

var_timeseries1 = []
for wrf_files in progressbar.progressbar(wrfoutfile_post):
    wrf_ncfile = Dataset(wrf_files)
    var = wrf_assign_coords(getvar(wrf_ncfile, var_name))
    var_timeseries1.append(
        var.interp(south_north=location[1],
                   west_east=location[0], method="nearest")
    )

var_timeseries_merged = xr.concat(var_timeseries, dim="time")
var_timeseries_merged1 = xr.concat(var_timeseries1, dim="time")

fig, axs = plt.subplots(1, 1, figsize=(7, 4))
axs.plot(var_timeseries_merged["Time"],
         var_timeseries_merged, 'b', label='LULC 2001')
axs.plot(var_timeseries_merged["Time"],
         var_timeseries_merged1, 'r', label='LULC 2017')
axs.set_xlabel("Time")
axs.set_ylabel(var_name)
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
