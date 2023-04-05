from netCDF4 import Dataset
import matplotlib.dates as mdates
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import wrf
import pandas as pd
from wrf import (
    getvar,
    interplevel,
    to_np,
    latlon_coords,
    get_cartopy,
    cartopy_xlim,
    cartopy_ylim,
    ALL_TIMES,
)
import glob
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import re
import os
import xarray as xr
import datetime
import progressbar


import tropycal.tracks as tracks

from scipy.interpolate import interp2d


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def convert_polar(var_cropped, margin=5, resolution=0.09):
    r = np.arange(0, margin, resolution)
    ang = np.arange(0, 361, 1) * np.pi / 180
    r_mesh, ang_mesh = np.meshgrid(r, ang)

    x_polar = r_mesh * np.cos(ang_mesh)
    y_polar = r_mesh * np.sin(ang_mesh)

    pcp_interp = interp2d(
        var_cropped["XLONG"].isel(south_north=1),
        var_cropped["XLAT"].isel(west_east=1),
        var_cropped.values,
    )

    pcp_polar = np.ones(x_polar.shape) * np.nan

    for x_id in range(y_polar.shape[0]):
        for y_id in range(y_polar.shape[1]):
            pcp_polar[x_id, y_id] = pcp_interp(
                x_polar[x_id, y_id], y_polar[x_id, y_id])

    return r, ang, pcp_polar


def wrf_assign_coords(var):
    lon_1d = var.XLONG.mean(dim="south_north").values
    lat_1d = var.XLAT.mean(dim="west_east").values
    return var.drop(["XLONG", "XLAT"]).assign_coords(
        south_north=lat_1d, west_east=lon_1d
    )


def get_precip(wrf_runs):
    wrfoutfile = sorted(glob.glob(wrf_runs + "wrfout_d01*-25_*"))

    basin = tracks.TrackDataset(
        basin="north_atlantic", source="hurdat", include_btk=False
    )

    harvey = basin.get_storm(("harvey", 2017))
    harvey_finer = harvey.to_xarray().interp(
        time=pd.date_range(
            harvey.to_xarray()["time"][0].values,
            harvey.to_xarray()["time"][-1].values,
            freq="1H",
        )
    )

    wrf_ncfile = Dataset(wrfoutfile[0])

    time_id = 30
    var_name = "slp"
    margin = 4
    ws_max = []
    rainfall = []
    wrf_time = []
    precip_polar = []
    wrf_ncfile_t0 = Dataset(wrfoutfile[0])

    for timeid in progressbar.progressbar(range(len(wrfoutfile) - 1)):

        wrf_ncfile2 = Dataset(wrfoutfile[timeid + 1])
        wrf_ncfile1 = Dataset(wrfoutfile[timeid])

        lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
        # var = getvar(wrf_ncfile, "uvmet10_wspd_wdir").sel(wspd_wdir="wspd")
        slp = wrf_assign_coords(getvar(wrf_ncfile2, "slp"))

        ref_var = getvar(wrf_ncfile, "uvmet10_wspd_wdir")
        wrf_time.append(ref_var["Time"].values)

        var_rainfall = (
            getvar(wrf_ncfile2, "RAINC") + getvar(wrf_ncfile2, "RAINNC")
        ) - (getvar(wrf_ncfile1, "RAINC") + getvar(wrf_ncfile1, "RAINNC"))

        #    var = getvar(wrf_ncfile, "RAINC")
        #    var.values = var_rainfall.values
        var = var_rainfall
        #    var = getvar(wrf_ncfile, "uvmet10_wspd_wdir").sel(wspd_wdir='wspd')*2

        cen_lon = harvey_finer.sel(time=ref_var["Time"].values)["lon"].values
        cen_lat = harvey_finer.sel(time=ref_var["Time"].values)["lat"].values

        cen_lon = slp.west_east[np.where(slp == slp.min())[1]].values
        cen_lat = slp.south_north[np.where(slp == slp.min())[0]].values

        lon_1d = ref_var.XLONG.mean(dim="south_north").values
        lat_1d = ref_var.XLAT.mean(dim="west_east").values

        lon_ids = np.where(
            np.logical_and(lon_1d >= cen_lon - margin,
                           lon_1d <= cen_lon + margin)
        )[0]

        lat_ids = np.where(
            np.logical_and(lat_1d >= cen_lat - margin,
                           lat_1d <= cen_lat + margin)
        )[0]

        var_cropped = var.isel(west_east=lon_ids, south_north=lat_ids)

        var_cropped = var_cropped.assign_coords(
            XLONG=var_cropped.XLONG - cen_lon)
        var_cropped = var_cropped.assign_coords(
            XLAT=var_cropped.XLAT - cen_lat)

        resolution = np.unique(np.round(np.diff(lon_1d), 2))[0]
        r, ang, var_polar = var_cropped_polar = convert_polar(var_cropped)
        precip_polar.append(var_polar)
    return r, ang, precip_polar


pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"


r, ang, pre_precip = get_precip(pre_files)
r, ang, post_precip = get_precip(post_files)

pre_precip_polar_mean = np.nansum(np.array(pre_precip), axis=0)
post_precip_polar_mean = np.nansum(np.array(post_precip), axis=0)

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='polar')
cont = ax.contourf(ang, r[:15]*111.11, post_precip_polar_mean[:, :15].T)
ax.set_theta_zero_location('N')
cbar = plt.colorbar(cont, ax=ax)
cbar.ax.set_ylabel('Precipitation (mm/hr)')
plt.title('LULC 2017')
plt.tight_layout()


fig = plt.figure(2)
ax = fig.add_subplot(111, projection='polar')
cont = ax.contourf(ang, r[:15]*111.11, pre_precip_polar_mean[:, :15].T)
ax.set_theta_zero_location('N')
cbar = plt.colorbar(cont, ax=ax)
cbar.ax.set_ylabel('Precipitation (mm/hr)')
plt.title('LULC 2001')

plt.tight_layout()


fig = plt.figure(3)
ax = fig.add_subplot(111, projection='polar')
cont = ax.contourf(ang, r[:15]*111.11, post_precip_polar_mean[:,
                   :15].T - pre_precip_polar_mean[:, :15].T, cmap='bwr')
ax.set_theta_zero_location('N')
cbar = plt.colorbar(cont, ax=ax)
cbar.ax.set_ylabel('Precipitation Error (mm/hr)')
plt.tight_layout()


plt.figure(figsize=(7, 4))
plt.plot(r * 111.111, np.nanmean(pre_precip_polar_mean.T, axis=1),
         'r', label='LULC 2001')
plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T, axis=1),
         'b', label='LULC 2017')
plt.xlabel("Radius from the storm center (km)")
plt.ylabel("Precipation")
plt.tight_layout()
plt.legend()
plt.xlim((0, np.round(r.max() * 111.11)))


plt.figure(figsize=(7, 4))
plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T -
         pre_precip_polar_mean.T, axis=1), 'r', label='LULC 2001')
# plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T, axis=1), 'b', label='LULC 2017')
plt.xlabel("Radius from the storm center (km)")
plt.ylabel("Precipation")
plt.tight_layout()
plt.legend()
plt.xlim((0, np.round(r.max() * 111.11)))
plt.show()
