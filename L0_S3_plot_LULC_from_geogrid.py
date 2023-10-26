import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import matplotlib.patches as patches
from osgeo import gdal, osr
import pyproj

from scipy.interpolate import griddata
from src.lulc_colormap import get_lulc_colormap
from src.coast import plot_coast


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

geog_file = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/pre/WPS/geo_em.d03.nc"
geog_file = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WPS/geo_em.d03.nc"
pre_geog_file = xr.open_dataset(geog_file, engine="netcdf4")
wrf_longitudes_pre = pre_geog_file["XLONG_M"].squeeze().values
wrf_latitudes_pre = pre_geog_file["XLAT_M"].squeeze().values
wrf_lulc_pre = pre_geog_file["LU_INDEX"].squeeze().squeeze().values


geog_file = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post/WPS/geo_em.d03.nc"
geog_file = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WPS/geo_em.d03.nc"
post_geog_file = xr.open_dataset(geog_file, engine="netcdf4")
wrf_longitudes_post = post_geog_file["XLONG_M"].squeeze().values
wrf_latitudes_post = post_geog_file["XLAT_M"].squeeze().values
wrf_lulc_post = post_geog_file["LU_INDEX"].squeeze().squeeze().values


def get_bb(file_name):
    return np.array(
        (
            file_name["XLONG_M"].min().values,
            file_name["XLONG_M"].max().values,
            file_name["XLAT_M"].min().values,
            file_name["XLAT_M"].max().values,
        )
    )


lulc_cmap, lulc_classes = get_lulc_colormap()


fig = plt.figure(figsize=(15, 7))
axes = plt.axes(projection=ccrs.PlateCarree())

img = axes.pcolormesh(
    wrf_longitudes_pre,
    wrf_latitudes_pre,
    wrf_lulc_pre,
    vmin=1,
    vmax=18,
    cmap=lulc_cmap,
    shading="auto",
)
plot_coast(axes)
cbar = plt.colorbar(img, ticks=np.arange(1, 18, 1) + 0.5, shrink=0.8)
cbar.ax.set_yticklabels(list(lulc_classes.keys()))  # horizontal colorbar
axes.set_extent(get_bb(pre_geog_file))
plt.title("LULC 2001")
plt.tight_layout()
plt.savefig('.'.join(geog_file.split("/")[-1].split(".")[:2])+'_LULC_pre_map.jpeg', dpi=300)

fig = plt.figure(figsize=(15, 7))
axes = plt.axes(projection=ccrs.PlateCarree())

img = axes.pcolormesh(
    wrf_longitudes_post,
    wrf_latitudes_post,
    wrf_lulc_post,
    vmin=1,
    vmax=18,
    cmap=lulc_cmap,
    shading="auto",
)
plot_coast(axes)
cbar = plt.colorbar(img, ticks=np.arange(1, 18, 1) + 0.5, shrink=0.8)
cbar.ax.set_yticklabels(list(lulc_classes.keys()))  # horizontal colorbar
axes.set_extent(get_bb(post_geog_file))
plt.title("LULC 2017")
plt.tight_layout()
plt.savefig('.'.join(geog_file.split("/")[-1].split(".")[:2])+'_LULC_post_map.jpeg', dpi=300)



def plot_lulc_bar(wrf_lulc_pre, wrf_lulc_post):
    assert wrf_lulc_pre.shape == wrf_lulc_post.shape == (17,)

    total_area_pre = np.sum(wrf_lulc_pre)
    total_area_post = np.sum(wrf_lulc_post)

    percentages_pre = (wrf_lulc_pre / total_area_pre) * 100
    percentages_post = (wrf_lulc_post / total_area_post) * 100

    lulc_classes = [
        "Class 1",
        "Class 2",
        "Class 3",
        "Class 4",
        "Class 5",
        "Class 6",
        "Class 7",
        "Class 8",
        "Class 9",
        "Class 10",
        "Class 11",
        "Class 12",
        "Class 13",
        "Class 14",
        "Class 15",
        "Class 16",
        "Class 17",
    ]

    width = 0.35
    x = np.arange(len(lulc_classes))

    fig, ax = plt.subplots(figsize=(14, 7))
    bars_pre = ax.bar(
        x - width / 2, percentages_pre, width, label="Pre Scenario", color="blue"
    )
    bars_post = ax.bar(
        x + width / 2, percentages_post, width, label="Post Scenario", color="red"
    )

    ax.set_xlabel("LULC Classes")
    ax.set_ylabel("Percentage")
    ax.set_title("Percentage of LULC Classes for Pre and Post Scenarios")
    ax.set_xticks(x)
    ax.set_xticklabels(lulc_classes, rotation=45, ha="right")
    ax.legend()

    def autolabel(bars):
        for bar in bars:
            height = bar.get_height()
            ax.annotate(
                f"{height:.2f}%",
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
            )

    autolabel(bars_pre)
    autolabel(bars_post)

    plt.tight_layout()


wrf_lulc_pre_count = np.array(
    [
        np.where(wrf_lulc_pre == lulc_classes_id)[0].shape[0]
        for lulc_classes_id in np.arange(1, 18)
    ]
)
wrf_lulc_post_count = np.array(
    [
        np.where(wrf_lulc_post == lulc_classes_id)[0].shape[0]
        for lulc_classes_id in np.arange(1, 18)
    ]
)

plot_lulc_bar(wrf_lulc_pre_count, wrf_lulc_post_count)

#plt.show()


pd.DataFrame(
    np.array((wrf_lulc_pre_count, wrf_lulc_post_count)).T,
    index=lulc_classes,
    columns=["pre", "post"],
).to_csv(
    "../tables/" + ".".join(geog_file.split("/")[-1].split(".")[:2]) + "_LULC_count.csv"
)



