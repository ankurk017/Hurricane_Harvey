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
import imageio

plt.rcParams.update({"font.size": 12, "font.weight": "bold"})

# matplotlib.use('Agg')

case = "_cntl"

urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 1.2}

urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 2.5}
domain_bb = [
    urban_change["west_east"] - urban_change["box"],
    urban_change["west_east"] + urban_change["box"],
    urban_change["south_north"] - urban_change["box"],
    urban_change["south_north"] + urban_change["box"],
]

levels = np.array(
    (0.1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700)
)
levels = np.array((0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200))
levels = np.array((0.01, 0.5, 1.0, 2, 4, 5, 7, 10, 15, 20, 40, 50, 70, 80, 90, 100))
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

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre/WRF{case}/test/em_real/wrfout_d02_2017-*")
)[
    37:52
]  # [7:12]# [22:24]
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post/WRF{case}/test/em_real/wrfout_d02_2017-*")
)[
    37:52
]  # [7:12] #[22:24]

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

prefiles = wrfoutfile_pre[:index][7:9]
postfiles = wrfoutfile_post[:index][7:9]

output_dir = f"../figures_paper/Ensemble/Rainfall/{case[1:]}/"


image_files = []
years = (
    2001,
    2017,
)
for fileid in progressbar.progressbar(range(len(prefiles) - 1)):
    wrf_coords = wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "RAINC"))
    wrf_pcp_pre = (
        getvar(Dataset(prefiles[fileid + 1]), "RAINC")
        + getvar(Dataset(prefiles[fileid + 1]), "RAINNC")
    ) - (
        getvar(Dataset(prefiles[fileid]), "RAINC")
        + getvar(Dataset(prefiles[fileid]), "RAINNC")
    )
    wrf_pcp_pre_uv = crop_region(
        wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "uvmet10")),
        location=urban_change,
    )

    # 	wrf_pcp_post = getvar(Dataset(postfiles[fileid+1]), 'wspd_wdir10').sel(wspd_wdir='wspd')
    wrf_pcp_post = (
        getvar(Dataset(postfiles[fileid + 1]), "RAINC")
        + getvar(Dataset(postfiles[fileid + 1]), "RAINNC")
    ) - (
        getvar(Dataset(postfiles[fileid]), "RAINC")
        + getvar(Dataset(postfiles[fileid]), "RAINNC")
    )
    wrf_pcp_post_uv = crop_region(
        wrf_assign_coords(getvar(Dataset(postfiles[fileid + 1]), "uvmet10")),
        location=urban_change,
    )

    slp = getvar(Dataset(prefiles[fileid + 1]), "slp")
    wrf_lat, wrf_lon = latlon_coords(slp)
    pre_post = (wrf_pcp_pre, wrf_pcp_post)
    pre_post_uv = (wrf_pcp_pre_uv, wrf_pcp_post_uv)

    fig, ax = plt.subplots(
        3,
        1,
        figsize=(4, 10),
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
        # 				      levels=levels, )
        ax[pre_post_id].quiver(
            pre_post_uv[pre_post_id]["west_east"][::10],
            pre_post_uv[pre_post_id]["south_north"][::10],
            pre_post_uv[pre_post_id].sel(u_v="u")[::10, ::10],
            pre_post_uv[pre_post_id].sel(u_v="v")[::10, ::10],
            units="x",
            scale=40,
            width=0.022,
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
        ax[pre_post_id].set_title(
            f"LULC {years[pre_post_id]} | " + f"{str(slp.Time.values)[:13]}", fontsize=12
        )
        cmap = copy.copy(cont.get_cmap())
        cmap.set_under("w")
        cont.set_cmap(cmap)

    #    cbar = plt.colorbar(cont, ax=ax[:2].ravel(), ticks=levels, shrink=0.75)
    #    cbar.ax.set_ylabel("Precipitation (mm/h)")

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
    ax[2].set_title(f"LULC 2017 - LULC 2001", fontsize=12)
    #    cbar = plt.colorbar(cont, ax=ax[2], shrink=0.75)
    #    cbar.ax.set_ylabel('Precipitation Error (mm/h)')
    plt.tight_layout()
    image_file = output_dir + f"/Houston_{str(slp.Time.values)[:13]}.jpeg"

plt.show()
