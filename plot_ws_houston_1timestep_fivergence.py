from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np
from scipy.interpolate import griddata
import progressbar
import matplotlib
from src.wrf_src import wrf_assign_coords, crop_region
from metpy.units import units
import metpy.calc as mpcalc

# matplotlib.use('Agg')

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


domain_bb = [-100, -93, 25.5, 31.5]
domain_bb = [-96.72, -94.21, 28.71, 30.91]
urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 1.2}
urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 0.5}

levels = np.array(
    (0.1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700)
)
levels = np.array((0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200))
levels = np.array(
    (
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.5,
        2,
        3,
        4,
        5,
        7,
        10,
        15,
        20,
        40,
        50,
    )
)
levels = np.array(
    (
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        1.5,
        2,
        3,
        4,
        5,
        7,
        10,
        15,
        20,
        40,
        50,
        70,
        100,
    )
)
levels = np.array((0.01, 0.5, 1.0, 2, 4, 5, 7, 10, 15, 20, 40, 50, 70, 80, 90, 100))
levels = np.array((1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150))
levels  = np.linspace(-4e-3, 4e-3, 25)
diff_levels = np.arange(-50, 55, 5)

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"
wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre_UCM/WRF/test/em_real/wrfout_d02_2017-*")
) #[22:24]
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF/test/em_real/wrfout_d02_2017-*")
) #[22:24]

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

prefiles = wrfoutfile_pre[:index]
postfiles = wrfoutfile_post[:index]


years = (
    2001,
    2017,
)
for fileid in progressbar.progressbar(range(len(prefiles) - 1)):
    # 	wrf_pcp_pre = getvar(Dataset(prefiles[fileid+1]), 'wspd_wdir10').sel(wspd_wdir='wspd')
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

    divergence_pre = mpcalc.divergence(wrf_pcp_pre_uv.sel(u_v="u").values, wrf_pcp_pre_uv.sel(u_v="v").values, dx=units.Quantity(0.027*111.11*1000, "m"), dy=units.Quantity(0.027*111.11*1000, "m"))
    divergence_post = mpcalc.divergence(wrf_pcp_post_uv.sel(u_v="u").values, wrf_pcp_post_uv.sel(u_v="v").values, dx=units.Quantity(0.027*111.11*1000, "m"), dy=units.Quantity(0.027*111.11*1000, "m"))

    slp = getvar(Dataset(prefiles[fileid + 1]), "slp")
    wrf_lat, wrf_lon = latlon_coords(slp)
    pre_post = (wrf_pcp_pre, wrf_pcp_post)
    pre_post_uv = (wrf_pcp_pre_uv, wrf_pcp_post_uv)
    pre_post_div = (divergence_pre, divergence_post)

    fig, ax = plt.subplots(   1,    3,     figsize=(18, 5),        sharey=True,        subplot_kw={"projection": ccrs.PlateCarree()},    )
    for pre_post_id in range(2):

        cont = ax[pre_post_id].contourf(
            pre_post_uv[pre_post_id]["west_east"],
            pre_post_uv[pre_post_id]["south_north"],
            pre_post_div[pre_post_id],
            cmap="PiYG",
            levels=levels,
       #     locator=ticker.LogLocator(),
        )
        # 				      levels=levels, )
        ax[pre_post_id].quiver(
            pre_post_uv[pre_post_id]["west_east"][::2],
            pre_post_uv[pre_post_id]["south_north"][::2],
            pre_post_uv[pre_post_id].sel(u_v="u")[::2, ::2] ,
            pre_post_uv[pre_post_id].sel(u_v="v")[::2, ::2],  units='x', scale=60, width=0.012, 
        )
        shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
        reader = shpreader.Reader(shapefile_path)
        geometries = reader.geometries()
        for geometry in geometries:
            ax[pre_post_id].add_geometries([geometry], ccrs.PlateCarree(),
        					 facecolor='none', edgecolor='blue', linewidth=0.5,)
       	coast.plot_coast(ax[pre_post_id])
        ax[pre_post_id].set_xlim((domain_bb[0], domain_bb[1]))
        ax[pre_post_id].set_ylim((domain_bb[2], domain_bb[3]))
        ax[pre_post_id].set_title(f"WRF (LULC {years[pre_post_id]})")
        # ax[1].set_yticklabels([])
        cbar = plt.colorbar(cont, ax=ax[pre_post_id], ticks=levels, shrink=0.75)
        cbar.ax.set_yticklabels(np.round(levels, 3))
        cbar.ax.set_ylabel("Divergence (1/s)")

    cont = ax[2].contourf(pre_post_uv[pre_post_id]["west_east"], pre_post_uv[pre_post_id]["south_north"], divergence_post-divergence_pre, cmap="seismic")
    shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()
    for geometry in geometries:
             ax[2].add_geometries([geometry], ccrs.PlateCarree(),
     				 facecolor='none', edgecolor='blue')
    coast.plot_coast(ax[2])
    ax[2].set_xlim((domain_bb[0], domain_bb[1]))
    ax[2].set_ylim((domain_bb[2], domain_bb[3]))
    ax[2].set_title(f'WRF (LULC {years[pre_post_id]})')
    cbar = plt.colorbar(cont, ax=ax[2], shrink=0.75)
    # 	cbar = plt.colorbar(cont, ax=ax[2], ticks=diff_levels, shrink=0.75)
    # 	cbar.ax.set_yticklabels(diff_levels)
    # 	cbar.ax.set_ylabel('Precipitation Error (mm)')
    plt.tight_layout()
    plt.suptitle(str(slp.Time.values)[:13], fontweight="bold")
# 	plt.show()
    plt.savefig(f"../figures/houston_pcp_div/Houston_{str(slp.Time.values)[:13]}.jpeg")
    plt.close()

#plt.show()





