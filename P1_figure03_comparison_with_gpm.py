from src.wrf_src import wrf_assign_coords
import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
#from self_utils import ahps, coast
import src.coast as coast
import numpy as np
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
import datetime
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/'

wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_cntl//test/em_real/wrfout_d03_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_cntl//test/em_real/wrfout_d03_2017-*'))


var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location

box = 1.5

var_timeseries_pre = []
var_timeseries_post = []
var_timeseries_pre_mean = []
var_timeseries_post_mean = []

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

    var_timeseries_pre.append(var_pre.sel(south_north=slice(location[1] - box, location[1] + box), west_east=slice(location[0] - box, location[0] + box)).mean())

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

    var_timeseries_post.append(var_post.sel(south_north=slice(location[1] - box, location[1] + box), west_east=slice(location[0] - box, location[0] + box)).mean())



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
gpm_rainfall_region = gpm.sel( lon=slice(location[0] - box, location[0] + box ),   lat=slice( location[1] - box, location[1] + box)         )

gpm_rainfall_region_mean = gpm_rainfall_region.mean(dim=["lon", "lat"])

#################

gpm_index = np.where([gpm_timestep_values==datetime.datetime.utcfromtimestamp(var_timeseries_post_merged["Time"][0].values.astype('O') / 1e9) for gpm_timestep_values in gpm_timestep])[0][0]
fig, axs = plt.subplots(1, 1, figsize=(7.5, 4.3))
axs.plot(gpm_timestep, np.cumsum(gpm_rainfall_region_mean.values)-np.cumsum(gpm_rainfall_region_mean.values)[gpm_index], "k-", label="GPM")
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
#plt.savefig('../figures/rainfall/WRF_GPM_rainfall_accum.jpeg')


fig, axs = plt.subplots(1, 1, figsize=(10, 5.5))
axs.plot(gpm_timestep, gpm_rainfall_region_mean.values, "k-", label="GPM")
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

plt.grid()
plt.minorticks_on()
axs.xaxis.set_minor_locator(MaxNLocator(46))

# plt.ylim((0, 510))
plt.tight_layout()
plt.savefig('../figures_paper/WRF_GPM_rainfall.jpeg')
plt.show()


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
#plt.savefig('../figures/rainfall/region_gpm_mean_cross.jpeg')
plt.show()



