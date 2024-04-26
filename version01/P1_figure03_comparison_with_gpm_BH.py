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
from matplotlib import ticker 


plt.rcParams.update({"font.size": 16, "font.weight": "bold"})

import numpy as np
import xarray as xr
from netCDF4 import Dataset
import progressbar  # You need to import the progressbar library

def get_wrf_pcp(wrfoutfile_pre_f, location_f, box=0.75):
    var_timeseries_pre_f = []

    for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_pre_f) - 1)):
        wrf_ncfile_pre = Dataset(wrfoutfile_pre_f[wrf_files_id + 1])
        wrf_ncfile_pre_tminus1 = Dataset(wrfoutfile_pre_f[wrf_files_id])

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
                var_pre["south_north"].values > location_f[1] - box,
                var_pre["south_north"].values <= location_f[1] + box,
            )
        )[0]
        lon_id = np.where(
            np.logical_and(
                var_pre["west_east"].values > location_f[0] - box,
                var_pre["west_east"].values <= location_f[0] + box,
            )
        )[0]

        var_timeseries_pre_f.append(var_pre.isel(south_north=lat_id, west_east=lon_id).mean())

    return xr.concat(var_timeseries_pre_f, dim="Time")


wrfoutfile_pbl01 = sorted(glob.glob('/rgroup/stela/akumar/BH_Simulations/2000/WRF_452_BH/test/em_real/wrfout_d02*'))
wrfoutfile_pbl02 = sorted(glob.glob('/rgroup/stela/akumar/BH_Simulations/2020/WRF_452_BH/test/em_real/wrfout_d02*'))
wrfoutfile_2050 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD//WRF_Building_Heights/WRF_Simulations/WRF_2050/WRF/test/em_real/wrfout_d02*'))
wrfoutfile_2100 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD//WRF_Building_Heights/WRF_Simulations/WRF_2100/WRF/test/em_real/wrfout_d02*'))



location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location
box=0.75


var_timeseries_pbl00 = get_wrf_pcp(wrfoutfile_pbl01, location, box=box)
var_timeseries_pbl01 = get_wrf_pcp(wrfoutfile_pbl02, location, box=box)
var_timeseries_pbl02 = get_wrf_pcp(wrfoutfile_2050, location, box=box)
var_timeseries_pbl03 = get_wrf_pcp(wrfoutfile_2100, location, box=box)

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
gpm_index = np.where([gpm_timestep_values==datetime.datetime.utcfromtimestamp(var_timeseries_pbl02["Time"][0].values.astype('O') / 1e9) for gpm_timestep_values in gpm_timestep])[0][0]
#################


fig, axs = plt.subplots(1, 1, figsize=(10, 5.5))
axs.plot(gpm_timestep, np.cumsum(gpm_rainfall_region.values)-np.cumsum(gpm_rainfall_region.values)[gpm_index], "k-", label="GPM")

axs.plot(var_timeseries_pbl00["Time"],  np.cumsum(var_timeseries_pbl00),    "b-",    label="BH 2001",)
axs.plot(var_timeseries_pbl01["Time"],  np.cumsum(var_timeseries_pbl01),    "r-",    label="BH 2020",)
axs.plot(var_timeseries_pbl02["Time"],  np.cumsum(var_timeseries_pbl02),    "g-",    label="BH 2050",)
axs.plot(var_timeseries_pbl03["Time"],  np.cumsum(var_timeseries_pbl03),    "m-",    label="BH 2100",)
#axs.plot(var_timeseries_pbl04["Time"],  np.cumsum(var_timeseries_pbl04),    "m-",    label="EDMF QNSE",)


# axs.plot([gpm_timestep[102], gpm_timestep[102]], [-10, 520], 'k--')
axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm)")
plt.legend()
plt.xticks(rotation=25)
plt.xlim(
    (var_timeseries_pbl00["Time"][0], var_timeseries_pbl00["Time"][-1])
)
axs.xaxis.set_minor_locator(ticker.MaxNLocator(46))
# plt.ylim((0, 510))
plt.tight_layout()
plt.savefig('../figures_paper/BH_projections/PBL_WRF_GPM_rainfall_accum.jpeg')
plt.show()

fig, axs = plt.subplots(1, 1, figsize=(10, 5.5))
axs.plot(gpm_timestep, gpm_rainfall_region.values, "k-", label="GPM")

axs.plot(var_timeseries_pbl00["Time"],  var_timeseries_pbl00,    "b-",    label="BH 2001",)
axs.plot(var_timeseries_pbl01["Time"],  var_timeseries_pbl01,    "r-",    label="BH 2017",)
axs.plot(var_timeseries_pbl02["Time"],  var_timeseries_pbl02,    "g-",    label="BH 2050",)
axs.plot(var_timeseries_pbl03["Time"],  var_timeseries_pbl02,    "m-",    label="BH 2100",)
#axs.plot(var_timeseries_pbl04["Time"],  var_timeseries_pbl04,    "m-",    label="EDMF QNSE",)

axs.set_xlabel("Time")
axs.set_ylabel("Precipitation (mm/hr)")
plt.legend()
plt.xticks(rotation=25)
plt.xlim(
    (var_timeseries_pbl00["Time"][0], var_timeseries_pbl00["Time"][-1])
)

plt.grid()
plt.minorticks_on()
axs.xaxis.set_minor_locator(ticker.MaxNLocator(46))

# plt.ylim((0, 510))
plt.tight_layout()
plt.savefig('../figures_paper/BH_projections/PBL_WRF_GPM_rainfall.jpeg')
plt.show()

start_lon, start_lat = location[0]-box, location[1]-box
end_lon, end_lat = location[0]+box, location[1]+box


fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
wrf_assign_coords(getvar(Dataset(wrfoutfile_pbl00[2]), "wspd_wdir10")).sel(wspd_wdir='wspd').plot(cmap='binary', levels=np.arange(0, 30, 2))
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
#plt.savefig('../figures_paper/PBL_WRF_GPM_rainfall_domain.jpeg')
plt.show()


