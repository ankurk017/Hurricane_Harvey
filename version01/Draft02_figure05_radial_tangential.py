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
    wa = getvar(ncfile, "wa", units="m/s")

    ua_pres = xr.concat([interplevel(ua, p, pres) for pres in pres_levels], dim="level")
    va_pres = xr.concat([interplevel(va, p, pres) for pres in pres_levels], dim="level")
    wa_pres = xr.concat([interplevel(wa, p, pres) for pres in pres_levels], dim="level")

    uv = xr.merge((ua_pres, va_pres, wa_pres))

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


def plot_rad_tang(uv_cs, xx="radius", case_name='pre', prefix=None, row = 0):

    # Radial Plot
    rad = axs[row, 0].contourf(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        cmap="bwr",
        levels=np.arange(-30, 40, 5),
    )
    contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(-60, 60, 10),
    )
    neg_contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(-20, 0, 5),
    )
    pos_contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, 40, 10),
    )
    axs[row, 0].clabel(contour_rad, inline=True, fmt="%1.1f", fontsize=12)

    skip = {'x': 4, 'y':4}

    rad_vals = uv_cs['wa_interp'] 
    xxx, yyy, = np.meshgrid(uv_cs[xx][::skip['x']], uv_cs['level'][::skip['y']])
    axs[row, 0].scatter(xxx, yyy, rad_vals[::skip['y'], ::skip['x']]*100, c=rad_vals[::skip['y'], ::skip['x']]*70, alpha=0.50, cmap='Greens', marker='o')
    axs[row, 0].scatter(xxx, yyy, -(rad_vals[::skip['y'], ::skip['x']]*500), c=rad_vals[::skip['y'], ::skip['x']]*70, alpha=1, cmap='Reds', marker='x')

    # Tangential Wind
    tan = axs[row, 1].contourf(uv_cs[xx], uv_cs["level"], uv_cs["tangential"], cmap="jet", levels=np.arange(-15, 40, 5))
    contour_tan = axs[row, 1].contour(
        uv_cs[xx], uv_cs["level"], uv_cs["tangential"], colors="k", linewidths=0.5
    )
    neg_contour_tan = axs[row, 1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(-20, 0, 5),
    )
    pos_contour_tan = axs[row, 1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, 50, 10),
    )
    axs[row, 1].clabel(contour_tan, inline=True, fmt="%1.1f", fontsize=12)
    #[ax.invert_yaxis() for ax in (axs[row, 0], axs[row, 1])]


    # Set color bar labels
    axs[row, 0].set_title(f"{case_name} | $V_r$ | {prefix}")
    axs[row, 1].set_title(f"{case_name} | $V_t$| {prefix}")

    [ax.set_xlabel("Distance from the storm center (km)") for ax in (axs[row, 0], axs[row, 1])]
    [ax.set_ylabel("Pressure (hPa)") for ax in (axs[row, 0], axs[row, 1])]

    if row == 0:
     cb_rad = plt.colorbar(rad, ax=axs[:, 0].ravel(), orientation='horizontal')
     cb_tan = plt.colorbar(tan, ax=axs[:, 1].ravel(), orientation='horizontal')
     cb_rad.set_label(f"Radial Wind ($V_r$)")
     cb_tan.set_label(f"Tangential Wind ($V_t$)")
     return (cb_rad, cb_tan)
    else:
     None

wrf_files_pre = sorted(
    glob.glob(
        f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_cntl/test/em_real/wrfout_d02_2017-08-27*"
    )
) #[6:13]
wrf_files_post = sorted(
    glob.glob(
        f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF_cntl/test/em_real/wrfout_d02_2017-08-27*"
    )
) #[6:13]


wrf_files_pre = sorted(glob.glob('/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-27*'))
wrf_files_post = sorted(glob.glob('/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-27*'))


outdir = '../figures_draft01/radial_tangential_FNL/'

for index in range(len(wrf_files_post)):
 print(index)
 prefix = wrf_files_post[index].split('/')[-1][11:24]

 rad_tang_pre = get_radial_tangential_wind(wrf_files_pre[index])
 rad_tang_post = get_radial_tangential_wind(wrf_files_post[index])


 fig, axs = plt.subplots(2, 2, figsize=(14, 9.5), sharex=True, sharey=True)

 rad_c, tan_c = plot_rad_tang(rad_tang_pre, case_name = 'LULC 2001', prefix = prefix, row = 0)
 plot_rad_tang(rad_tang_post, case_name = 'LULC 2017', prefix = prefix, row = 1)
 plt.subplots_adjust(left=0.09, right=0.95, wspace=0.12, hspace=0.25, bottom=0.17, top=0.9)

 rad_c.ax.set_position([0.15, 0.07, 0.3, 0.02])  
 tan_c.ax.set_position([0.58, 0.07, 0.3, 0.02])  
 axs[0, 1].invert_yaxis()

 plt.savefig(f"{outdir}/{prefix}.jpeg", dpi=400)
 plt.close()







