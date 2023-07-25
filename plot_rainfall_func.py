from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import src.coast as coast
import numpy as np
from scipy.interpolate import griddata
import progressbar
import matplotlib
from src.wrf_src import wrf_assign_coords, crop_region
from metpy.units import units
import metpy.calc as mpcalc
import matplotlib.colors as mcolors
import copy

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def plot_rainfall_and_winds(
    pre=None,
    post=None,
    start_point=None,
    end_point=None,
    plot_cross=False,
    case="_cntl",
):
    domain_bb = [-96.72, -94.21, 28.71, 30.91]
    urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 0.8}
    domain_bb = [
        urban_change["west_east"] - urban_change["box"],
        urban_change["west_east"] + urban_change["box"],
        urban_change["south_north"] - urban_change["box"],
        urban_change["south_north"] + urban_change["box"],
    ]

    levels = np.array(
        (
            2,
            10,
            20,
            30,
            40,
            50,
            60,
            70,
            80,
            90,
            100,
        )
    )
    diff_levels = np.arange(-50, 55, 5)

    years = (
        2001,
        2017,
    )
    for fileid in progressbar.progressbar(range(len(pre) - 1)):
        wrf_coords = wrf_assign_coords(getvar(Dataset(pre[fileid + 1]), "RAINC"))
        wrf_pcp_pre = (
            getvar(Dataset(pre[fileid + 1]), "RAINC")
            + getvar(Dataset(pre[fileid + 1]), "RAINNC")
        ) - (
            getvar(Dataset(pre[fileid]), "RAINC")
            + getvar(Dataset(pre[fileid]), "RAINNC")
        )
        wrf_pcp_pre_uv = crop_region(
            wrf_assign_coords(getvar(Dataset(pre[fileid + 1]), "uvmet10")),
            location=urban_change,
        )

        #  wrf_pcp_post = getvar(Dataset(postfiles[fileid+1]), 'wspd_wdir10').sel(wspd_wdir='wspd')
        wrf_pcp_post = (
            getvar(Dataset(post[fileid + 1]), "RAINC")
            + getvar(Dataset(post[fileid + 1]), "RAINNC")
        ) - (
            getvar(Dataset(post[fileid]), "RAINC")
            + getvar(Dataset(post[fileid]), "RAINNC")
        )
        wrf_pcp_post_uv = crop_region(
            wrf_assign_coords(getvar(Dataset(post[fileid + 1]), "uvmet10")),
            location=urban_change,
        )

        slp = getvar(Dataset(pre[fileid + 1]), "slp")
        wrf_lat, wrf_lon = latlon_coords(slp)
        pre_post = (wrf_pcp_pre, wrf_pcp_post)
        pre_post_uv = (wrf_pcp_pre_uv, wrf_pcp_post_uv)

        fig, ax = plt.subplots(
            1,
            3,
            figsize=(18, 5),
            sharey=True,
            subplot_kw={"projection": ccrs.PlateCarree()},
        )
        for pre_post_id in range(2):

            cont = ax[pre_post_id].contourf(
                wrf_coords["west_east"],
                wrf_coords["south_north"],
                pre_post[pre_post_id],
                cmap="YlOrRd",
                levels=levels,
                extend="both"
                #     locator=ticker.LogLocator(),
            )
            #              levels=levels, )
            ax[pre_post_id].quiver(
                pre_post_uv[pre_post_id]["west_east"][::4],
                pre_post_uv[pre_post_id]["south_north"][::4],
                pre_post_uv[pre_post_id].sel(u_v="u")[::4, ::4],
                pre_post_uv[pre_post_id].sel(u_v="v")[::4, ::4],
                units="x",
                scale=60,
                width=0.012,
            )
            shapefile_path = (
                "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
            )
            reader = shpreader.Reader(shapefile_path)
            geometries = reader.geometries()
            for geometry in geometries:
                ax[pre_post_id].add_geometries(
                    [geometry],
                    ccrs.PlateCarree(),
                    facecolor="none",
                    edgecolor="blue",
                    linewidth=0.5,
                )
            coast.plot_coast(ax[pre_post_id])
            ax[pre_post_id].set_xlim((domain_bb[0], domain_bb[1]))
            ax[pre_post_id].set_ylim((domain_bb[2], domain_bb[3]))
            ax[pre_post_id].set_title(f"WRF (LULC {years[pre_post_id]})")
            # ax[1].set_yticklabels([])
            cmap = copy.copy(cont.get_cmap())
            cmap.set_under("w")
            cont.set_cmap(cmap)

            cbar = plt.colorbar(cont, ax=ax[pre_post_id], ticks=levels, shrink=0.75)
            #        cbar.ax.set_yticklabels(np.round(levels, 3))
            cbar.ax.set_ylabel("Precipitation (mm/h)")

        cont = ax[2].contourf(
            wrf_coords["west_east"],
            wrf_coords["south_north"],
            pre_post[1] - pre_post[0],
            levels=diff_levels,
            cmap="bwr",
            extend="both",
        )
        shapefile_path = (
            "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
        )
        reader = shpreader.Reader(shapefile_path)
        geometries = reader.geometries()
        for geometry in geometries:
            ax[2].add_geometries(
                [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
            )
        coast.plot_coast(ax[2])
        ax[2].set_xlim((domain_bb[0], domain_bb[1]))
        ax[2].set_ylim((domain_bb[2], domain_bb[3]))
        ax[2].set_title(f"WRF (LULC {years[pre_post_id]})")
        cbar = plt.colorbar(cont, ax=ax[2], shrink=0.75)
        #  cbar = plt.colorbar(cont, ax=ax[2], ticks=diff_levels, shrink=0.75)
        #  cbar.ax.set_yticklabels(diff_levels)
        cbar.ax.set_ylabel("Precipitation Error (mm/h)")
        plt.tight_layout()
        plt.suptitle(
            f"{case[1:].upper()} | " + str(slp.Time.values)[:13], fontweight="bold"
        )
        if plot_cross:
            [ax.plot((end_point.lon, start_point.lon), (end_point.lat, start_point.lat), "g-", linewidth=2) for ax in ax]
    return None
