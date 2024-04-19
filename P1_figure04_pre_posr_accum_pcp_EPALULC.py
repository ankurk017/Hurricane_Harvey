from src.wrf_src import wrf_assign_coords
import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import src.coast as coast
import numpy as np
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

import datetime
from matplotlib.patches import Rectangle

plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"

wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre/WRF/test/em_real/wrfout_d02_2017-*")
)
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post/WRF//test/em_real/wrfout_d02_2017-*")
)

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))
index = 85
wrfoutfile_pre = wrfoutfile_pre[36 : 36 + 48]
wrfoutfile_post = wrfoutfile_post[36 : 36 + 48]


# comment this for multiple days

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations/"

dates = ["26", "27", "28"]
postfix = "2100"


wrfoutfile_pre, wrfoutfile_post = [], []
for date in dates:
    wrfoutfile_pre += sorted(
        glob.glob(home + f"/WRF_FNL_pre/WRF/test/em_real/wrfout_d02_2017-08-{date}*")
    )
    wrfoutfile_post += sorted(
        glob.glob(home + f"/WRF_FNL_{postfix}/WRF/test/em_real/wrfout_d02_2017-08-{date}*")
    )


var_name = "slp"

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499, 29.74)  # Houston location

box = 0.45  # works fine
box = 0.75
box = 1.25
box = 2.25

var_timeseries_pre = []
var_timeseries_post = []

# for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_post) - 1)):
for wrf_files_id in range(1):
    wrf_ncfile_pre = Dataset(wrfoutfile_pre[-1])
    wrf_ncfile_pre_tminus1 = Dataset(wrfoutfile_pre[0])

    var_pre_tmp = (
        wrf_assign_coords(
            getvar(wrf_ncfile_pre, "RAINC") + getvar(wrf_ncfile_pre, "RAINNC")
        ).values
        - wrf_assign_coords(
            getvar(wrf_ncfile_pre_tminus1, "RAINC")
            + getvar(wrf_ncfile_pre_tminus1, "RAINNC")
        ).values
    )

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

    var_timeseries_pre = var_pre.isel(south_north=lat_id, west_east=lon_id)

    wrf_ncfile_post = Dataset(wrfoutfile_post[-1])
    wrf_ncfile_post_tminus1 = Dataset(wrfoutfile_post[0])
    var_post_tmp = (
        wrf_assign_coords(
            getvar(wrf_ncfile_post, "RAINC") + getvar(wrf_ncfile_post, "RAINNC")
        ).values
        - wrf_assign_coords(
            getvar(wrf_ncfile_post_tminus1, "RAINC")
            + getvar(wrf_ncfile_post_tminus1, "RAINNC")
        ).values
    )

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

    var_timeseries_post = var_post.isel(south_north=lat_id, west_east=lon_id)


start_lon, start_lat = location[0] - box, location[1] - box
end_lon, end_lat = location[0] + box, location[1] + box


fig, axes = plt.subplots(
    nrows=1, ncols=3, figsize=(18, 5), subplot_kw={"projection": ccrs.PlateCarree()}
)

ax = axes[0]
var_timeseries_pre.plot(
    ax=ax,
    cmap="gist_ncar",
    levels=np.arange(0, 1050, 50),
    extend="both",
    cbar_kwargs={"shrink": 0.8},
)
ax.set_xlim([start_lon, end_lon])
rect = Rectangle(
    (start_lon, start_lat),
    end_lon - start_lon,
    end_lat - start_lat,
    facecolor="none",
    edgecolor="red",
    linewidth=2,
)
ax.add_patch(rect)
coast.plot_coast(ax)
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()
for geometry in geometries:
    ax.add_geometries(
        [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="black"
    )

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([start_lon, end_lon])
ax.set_ylim([start_lat, end_lat])
ax.set_title(f"LULC 2020")

ax = axes[1]
var_timeseries_post.plot(
    ax=ax,
    cmap="gist_ncar",
    levels=np.arange(0, 1050, 50),
    extend="both",
    cbar_kwargs={"shrink": 0.8},
)
ax.set_xlim([start_lon, end_lon])
rect = Rectangle(
    (start_lon, start_lat),
    end_lon - start_lon,
    end_lat - start_lat,
    facecolor="none",
    edgecolor="red",
    linewidth=2,
)
ax.add_patch(rect)
coast.plot_coast(ax)
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()
for geometry in geometries:
    ax.add_geometries(
        [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="black"
    )

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([start_lon, end_lon])
ax.set_ylim([start_lat, end_lat])
ax.set_title("LULC 2100")

ax = axes[2]
(var_timeseries_post - var_timeseries_pre).plot(
    ax=ax,
    cmap="bwr",
    levels=np.arange(-120, 140, 20),
    extend="both",
    cbar_kwargs={"shrink": 0.8},
)
ax.set_xlim([start_lon, end_lon])
rect = Rectangle(
    (start_lon, start_lat),
    end_lon - start_lon,
    end_lat - start_lat,
    facecolor="none",
    edgecolor="red",
    linewidth=2,
)
ax.add_patch(rect)
coast.plot_coast(ax)
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()
for geometry in geometries:
    ax.add_geometries(
        [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="black"
    )

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([start_lon, end_lon])
ax.set_ylim([start_lat, end_lat])

ax.set_title("Difference")

plt.tight_layout()
plt.savefig(
    f"../figures_paper/rainfall_Houston_all_accumulated_difference_{date}_EPALULC_{postfix}.jpeg"
)
# plt.show()


plt.rcParams.update({"font.size": 10, "font.weight": "bold"})
data_post = var_timeseries_post.values.ravel()
data_pre = var_timeseries_pre.values.ravel()

fig, ax = plt.subplots(3, 1, figsize=(5, 3.1), sharex=True)
ax[1].boxplot(
    data_post,
    vert=False,
    positions=[1],
    widths=0.6,
    showfliers=False,
    labels=["LULC 2017"],
)
ax[1].scatter(
    data_post.mean(), [1], marker="D", color="red", s=60, label="Mean LULC 2017"
)
# ax[0].set_xlabel('WRF accumulated precipitation (mm/hr)')


ax[0].boxplot(
    data_pre,
    vert=False,
    positions=[1],
    widths=0.6,
    showfliers=False,
    labels=["LULC 2001"],
)
ax[0].scatter(
    data_pre.mean(), [1], marker="D", color="red", s=60, label="Mean LULC 2001"
)
# ax[1].set_xlabel('WRF accumulated precipitation (mm/hr)')

data_diff = data_post - data_pre
ax[2].boxplot(
    data_diff,
    vert=False,
    positions=[1],
    widths=0.6,
    showfliers=False,
    labels=["2017 - 2001"],
)
ax[2].scatter(
    data_diff.mean(), [1], marker="D", color="red", s=60, label="Mean LULC 2001"
)
ax[2].set_xlabel("WRF accumulated precipitation error (mm/hr)")
plt.tight_layout()


plt.savefig(
    f"../figures_paper/rainfall_Houston_pre_post_boxplot_{date}_EPALULC.jpeg", dpi=300
)
# plt.show()
plt.close()
