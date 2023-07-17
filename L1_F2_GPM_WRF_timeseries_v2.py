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
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from self_utils import coast

import datetime
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


#wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
#wrfoutfile_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

#wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
#wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'
#wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/WRF_mov1_GFS_IC25_12UTC_v2/pre/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_2017-*'))
#wrfoutfile_post = sorted(glob.glob(home_2512 + f'/WRF_mov1_GFS_IC25_12UTC_v2/post/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_2017-*'))


wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/WRF_FNL_2612/WRF/test/em_real/wrfout_d01_2017-*')) # gives good results
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/WRF_FNL_2612/WRF/test/em_real/wrfout_d02_2017-*')) # gives good results



home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_Simulations/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/WRF_HRRR_sim01_pre/WRF/test/em_real/wrfout_d01_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/WRF_HRRR_sim01_post/WRF/test/em_real/wrfout_d01_2017-*'))

#wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/WRFV45_FNL_2518/WRF_old/test/em_real//wrfout_d01_2017-*'))
#wrfoutfile_post = wrfoutfile_pre

#home_nomoving = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'
#wrfoutfile_pre = sorted(glob.glob(home_nomoving + f'WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06/wrfout_d01*'))
#wrfoutfile_post = sorted(glob.glob(home_nomoving + f'WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06/wrfout_d01*'))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_1dom/test/em_real/wrfout_d01_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_1dom/test/em_real/wrfout_d01_2017-*'))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_2dom/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))

wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre_UCM/WRF//test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post_UCM/WRF//test/em_real/wrfout_d02_2017-*'))



var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location

box = 0.6

var_timeseries_pre = []
var_timeseries_post = []

for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_post) - 1)):

    wrf_ncfile_pre = Dataset(wrfoutfile_pre[wrf_files_id + 1])
    wrf_ncfile_pre_tminus1 = Dataset(wrfoutfile_pre[wrf_files_id])

    var_pre_tmp = wrf_assign_coords(
        getvar(wrf_ncfile_pre, "RAINC") + getvar(wrf_ncfile_pre, "RAINNC")
    ).values - wrf_assign_coords(
        getvar(wrf_ncfile_pre_tminus1, "RAINC")
        + getvar(wrf_ncfile_pre_tminus1, "RAINNC")
    ).values

    var_pre = wrf_assign_coords(getvar(wrf_ncfile_pre, "RAINC"))
    var_pre.values = var_pre_tmp

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
    ).values - wrf_assign_coords(
        getvar(wrf_ncfile_post_tminus1, "RAINC")
        + getvar(wrf_ncfile_post_tminus1, "RAINNC")
    ).values

    var_post = wrf_assign_coords(getvar(wrf_ncfile_post, "RAINC"))
    var_post.values = var_post_tmp
    lat_id = np.where(
        np.logical_and(
            var_post["south_north"].values > location[1] - box,
            var_post["south_north"].values <= location[1] + box,
        )
    )[0]
    lon_id = np.where(
        np.logical_and(
            var_post["west_east"].values > location[0] - box,
            var_post["west_east"].values <= location[0] + box,
        )
    )[0]


    var_timeseries_post.append(
        var_post.isel(south_north=lat_id, west_east=lon_id).mean()
    )

var_timeseries_pre_merged = xr.concat(var_timeseries_pre, dim="Time")
var_timeseries_post_merged = xr.concat(var_timeseries_post, dim="Time")

####### GPM ####
gpm_files = sorted(
    glob.glob("/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/GPM/*nc4")
#)[432:650]
)[432:750]
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

gpm_index = np.where([gpm_timestep_values==datetime.datetime.utcfromtimestamp(var_timeseries_post_merged["Time"][0].values.astype('O') / 1e9) for gpm_timestep_values in gpm_timestep])[0][0]
fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(gpm_timestep, np.cumsum(gpm_rainfall_region.values)-np.cumsum(gpm_rainfall_region.values)[gpm_index], "k-", label="GPM")
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
plt.savefig('../figures/rainfall/WRF_GPM_rainfall_accum.jpeg')


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
plt.savefig('../figures/rainfall/WRF_GPM_rainfall.jpeg')

start_lon, start_lat = location[0]-box, location[1]-box
end_lon, end_lat = location[0]+box, location[1]+box


fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
wrf_assign_coords(getvar(Dataset(wrfoutfile_pre[2]), "wspd_wdir10")).sel(wspd_wdir='wspd').plot(cmap='binary', levels=np.arange(0, 30, 2))
rect = Rectangle((start_lon, start_lat), end_lon - start_lon, end_lat - start_lat,
                 facecolor='none', edgecolor='red', linewidth=2)
ax.add_patch(rect)
coast.plot_coast(ax)
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()
for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([-99, -92.5])
ax.set_ylim([27, 32.7])

plt.tight_layout()
plt.savefig('../figures/rainfall/region_gpm_mean_cross.jpeg')
plt.show()



