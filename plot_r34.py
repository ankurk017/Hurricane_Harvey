from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from wrf import (
    getvar,
    latlon_coords,
)
import glob
import progressbar

from scipy.interpolate import interp2d

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def convert_polar(var_cropped, margin=5, resolution=0.09):
    r = np.arange(0, margin, resolution)
    ang = np.arange(0, 361, 1) * np.pi / 180
    r_mesh, ang_mesh = np.meshgrid(r, ang)

    x_polar = r_mesh * np.cos(ang_mesh)
    y_polar = r_mesh * np.sin(ang_mesh)

    pcp_interp = interp2d(
        var_cropped.west_east,
        var_cropped.south_north,
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


pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"


wrfoutfile = sorted(glob.glob(pre_files + "wrfout_d02*-25_*"))

time_id = 30
var_name = "slp"
margin = 3
r34 = []
wrf_time = []
precip_polar = []
for timeid in progressbar.progressbar(range(len(wrfoutfile) - 1)):

    wrf_ncfile2 = Dataset(wrfoutfile[timeid + 1])

    lats, lons = latlon_coords(getvar(wrf_ncfile2, "RAINC"))
    slp = wrf_assign_coords(getvar(wrf_ncfile2, "slp"))

    ref_var = getvar(wrf_ncfile2, "uvmet10_wspd_wdir", units="knots")
    wrf_time.append(ref_var["Time"].values)

    var = wrf_assign_coords(ref_var.sel(wspd_wdir="wspd"))

    cen_lon = var.west_east[np.where(var == var.min())[0]].values
    cen_lat = var.south_north[np.where(var == var.min())[1]].values

    lat_1d = var.south_north.values
    lon_1d = var.west_east.values

    lon_ids = np.where(
        np.logical_and(lon_1d >= cen_lon - margin, lon_1d <= cen_lon + margin)
    )[0]

    lat_ids = np.where(
        np.logical_and(lat_1d >= cen_lat - margin, lat_1d <= cen_lat + margin)
    )[0]

    var_cropped = var.isel(west_east=lon_ids, south_north=lat_ids)

    var_cropped = var_cropped.assign_coords(
        west_east=var_cropped.west_east - cen_lon)
    var_cropped = var_cropped.assign_coords(
        south_north=var_cropped.south_north - cen_lat
    )

    resolution = np.unique(np.round(np.diff(lon_1d), 2))[0]
    r, ang, var_polar = var_cropped_polar = convert_polar(
        var_cropped, margin=margin)
    print(" \n")
    precip_polar.append(var_polar)
    r34.append(
        np.nanmean(
            np.array(
                [
                    r[np.argmin(np.abs(np.array(var_polar[angles, :]) - 34))]
                    for angles in range(len(ang))
                ]
            )
        )
        * 111.11
    )


r34_values = pd.DataFrame(
    np.array([[str(wrf_time)[:13] for wrf_time in wrf_time], r34]).T)
print(r34_values)
