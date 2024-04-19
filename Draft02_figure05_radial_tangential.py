from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

from wrf import (
    to_np,
    getvar,
    smooth2d,
    get_cartopy,
    cartopy_xlim,
    cartopy_ylim,
    latlon_coords,
    interplevel,
)

# from self_utils import ahps, coast
import src.coast as coast
import src.ahps as ahps
import numpy as np
from scipy.interpolate import griddata
from matplotlib import gridspec
from src.wrf_src import wrf_assign_coords
import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
import xarray as xr
import metpy
from metpy.interpolate import cross_section
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import progressbar

plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "font.family": "monospace"}
)

def get_radial_tangential_wind(
    infile, pres_levels=np.arange(1000, 100, -10), radial_extent=5
):

    ncfile = Dataset(infile)

    slp = wrf_assign_coords(getvar(ncfile, "slp"))
    cen_lat = slp.south_north[np.where(slp == slp.min())[0]].values
    cen_lon = slp.west_east[np.where(slp == slp.min())[1]].values
    cen_loc = (cen_lon, cen_lat)

    slp = getvar(ncfile, "slp")
    p = getvar(ncfile, "pressure")

    ua = getvar(ncfile, "ua", units="m/s")
    va = getvar(ncfile, "va", units="m/s")

    ua_pres = xr.concat([interplevel(ua, p, pres) for pres in pres_levels], dim="level")
    va_pres = xr.concat([interplevel(va, p, pres) for pres in pres_levels], dim="level")
    uv = xr.merge((ua_pres, va_pres))

    uv_xr = wrf_assign_coords(
        uv.metpy.assign_crs(
            grid_mapping_name="latitude_longitude", earth_radius=6371229.0
        )
    )
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})

    uv_xr["y"] = uv_xr["y"].values - cen_loc[1]
    uv_xr["x"] = uv_xr["x"].values - cen_loc[0]

    start = (0, 0)
    theta = np.linspace(0, np.pi / 4, 8)
    # theta = np.linspace(0, 2*np.pi, 16)

    x = radial_extent * np.cos(theta)
    y = radial_extent * np.sin(theta)

    uv_rad_tang = []
    for index in progressbar.progressbar(range(len(x))):
        uv_cs = cross_section(uv_xr, start, (y[index], x[index])).set_coords(("y", "x"))
        uv_cs["radial"], uv_cs["tangential"] = metpy.calc.cross_section_components(
            uv_cs["ua_interp"], uv_cs["va_interp"]
        )
        uv_cs["radius"] = (np.sqrt(uv_cs["y"] ** 2 + uv_cs["x"] ** 2)) * 111.11
        uv_rad_tang.append(uv_cs)

    uv_radial_tangential_wind = xr.concat(uv_rad_tang, dim="new").mean(dim="new")
    return uv_radial_tangential_wind


def plot_rad_tang(uv_cs, xx="radius"):

    fig, axs = plt.subplots(1, 2, figsize=(16, 6.5), )

    # Radial Plot
    rad = axs[0].contourf(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        cmap="bwr",
        levels=np.arange(-30, 40, 5),
    )
    contour_rad = axs[0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(-60, 60, 10),
    )
    neg_contour_rad = axs[0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(uv_cs["radial"].min(), 0, 5),
    )
    pos_contour_rad = axs[0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, uv_cs["radial"].max(), 10),
    )
    axs[0].clabel(contour_rad, inline=True, fmt="%1.1f", fontsize=12)
    cb_rad = plt.colorbar(rad, ax=axs[0])

    # Tangential Wind
    tan = axs[1].contourf(uv_cs[xx], uv_cs["level"], uv_cs["tangential"], cmap="jet")
    contour_tan = axs[1].contour(
        uv_cs[xx], uv_cs["level"], uv_cs["tangential"], colors="k", linewidths=0.5
    )
    neg_contour_tan = axs[1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(uv_cs["tangential"].min(), 0, 5),
    )
    pos_contour_tan = axs[1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, uv_cs["tangential"].max(), 10),
    )
    axs[1].clabel(contour_tan, inline=True, fmt="%1.1f", fontsize=12)
    cb_tan = plt.colorbar(tan, ax=axs[1])

    [ax.invert_yaxis() for ax in axs]

    # Set color bar labels
    cb_rad.set_label(f"Radial Wind ($V_r$)")
    cb_tan.set_label(f"Tangential Wind ($V_t$)")

    # Set color bar labels
    axs[0].set_title(f"Radial Wind ($V_r$)")
    axs[1].set_title(f"Tangential Wind ($V_t$)")

    [ax.set_xlabel("Distance from the storm center (km)") for ax in axs]
    [ax.set_ylabel("Pressure (hPa)") for ax in axs]

    plt.tight_layout()
    # plt.savefig('test.jpeg')
    #plt.show()

case_name = 'post'
wrf_files = sorted(
    glob.glob(
        f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/{case_name}/WRF_cntl/test/em_real/wrfout_d02_2017-08-27*"
    )
)
outdir = '../figures_paper/radial_tangential/'

for input_file in wrf_files:
 print(input_file)
 rad_tang = get_radial_tangential_wind(input_file)
 
 plot_rad_tang(rad_tang)
 plt.suptitle(f"{case_name} | {input_file.split('/')[-1][:24]}")
 plt.tight_layout()
 figname = f"{outdir}/{case_name}/{case_name}_{input_file.split('/')[-1][:24]}.jpeg"
 plt.savefig(figname, dpi=400)
 plt.close()

#plt.show()







