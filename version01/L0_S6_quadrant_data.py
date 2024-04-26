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


wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_pre/WRF_9-3-1-51_levels/WRFV4/wrfout/"

wrfoutfile = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

basin = tracks.TrackDataset(basin="north_atlantic",
                            source="hurdat", include_btk=False)

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
wrf_ncfile_t0 = Dataset(wrfoutfile[0])

for timeid in progressbar.progressbar(range(len(wrfoutfile))):
    wrf_ncfile = Dataset(wrfoutfile[timeid])

    lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
    # var = getvar(wrf_ncfile, "uvmet10_wspd_wdir").sel(wspd_wdir="wspd")
    ref_var = getvar(wrf_ncfile, "uvmet10_wspd_wdir")
    wrf_time.append(ref_var["Time"].values)
    var_rainfall = (getvar(wrf_ncfile, "RAINC") + getvar(wrf_ncfile, "RAINNC")) - (
        getvar(wrf_ncfile_t0, "RAINC") + getvar(wrf_ncfile_t0, "RAINNC")
    )

    var = getvar(wrf_ncfile, "RAINC") + getvar(wrf_ncfile, "RAINNC")
    var.values = var_rainfall.values

    #    var = getvar(wrf_ncfile, "uvmet10_wspd_wdir").sel(wspd_wdir='wspd')*2

    cen_lon = harvey_finer.sel(time=ref_var["Time"].values)["lon"].values
    cen_lat = harvey_finer.sel(time=ref_var["Time"].values)["lat"].values

    lon_1d = ref_var.XLONG.mean(dim="south_north").values
    lat_1d = ref_var.XLAT.mean(dim="west_east").values

    lon_ids = np.where(
        np.logical_and(lon_1d >= cen_lon - margin, lon_1d <= cen_lon + margin)
    )[0]

    lat_ids = np.where(
        np.logical_and(lat_1d >= cen_lat - margin, lat_1d <= cen_lat + margin)
    )[0]

    var_cropped = var.isel(west_east=lon_ids, south_north=lat_ids)

    var_cropped = var_cropped.assign_coords(XLONG=var_cropped.XLONG - cen_lon)
    var_cropped = var_cropped.assign_coords(XLAT=var_cropped.XLAT - cen_lat)

    resolution = np.unique(np.round(np.diff(lon_1d), 2))[0]
    r, ang, var_polar = var_cropped_polar = convert_polar(var_cropped)
    out = []
    for angles in np.array((0, 90, 180, 270)):
        out.append(
            np.round(
                np.array(
                    [
                        np.nanmean(
                            var_polar[:, np.logical_and(r >= rad1, r <= rad2)].mean(-1)[
                                np.logical_and(
                                    ang * 180 / np.pi >= angles,
                                    ang * 180 / np.pi <= angles + 90,
                                )
                            ]
                        )
                        for rad1, rad2 in zip(np.array((0, 1, 2)), np.array((1, 2, 4)))
                    ]
                ),
                3,
            )
        )
    rainfall.append(np.array(out).T.ravel())

rainfall_quadrant = pd.DataFrame(
    np.array(rainfall),
    index=wrf_time,
    columns=[
        "RF (0-100 km)",
        "RR (0-100 km)",
        "LR (0-100 km)",
        "LF (0-100 km)",
        "RF (100-200 km)",
        "RR (100-200 km)",
        "LR (100-200 km)",
        "LF (100-200 km)",
        "RF (200-400 km)",
        "RR (200-400 km)",
        "LR (200-400 km)",
        "LF (200-400 km)",
    ],
)
rainfall_quadrant.to_csv('Rainfall_Harvey_quadrant.csv')

"""
plt.figure()
fig = plt.axes(projection='polar')
plt.contourf(ang, r, pcp_polar.T)
plt.show()
plt.plot(r*111.111, np.nanmean(pcp_polar.T, axis=1))
plt.xlabel('Radius from the storm center (km)')
plt.ylabel('Precipation')
plt.show()
"""
