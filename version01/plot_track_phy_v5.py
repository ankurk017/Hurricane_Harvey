from self_utils import coast
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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import re
import os
import xarray as xr
import datetime
import progressbar

import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep

import tropycal.tracks as tracks

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def get_track_details(wrf_files):
    track_lon = []
    track_lat = []
    slp_min = []
    ws_max = []
    for timeid in progressbar.progressbar(range(len(wrf_files))):
        wrf_ncfile = Dataset(wrf_files[timeid])
        slp = getvar(wrf_ncfile, "slp")
        slp_min.append(slp.min())
        ws = getvar(wrf_ncfile, "uvmet10_wspd_wdir").sel(wspd_wdir="wspd")
        ws_max.append(ws.max())
        lats, lons = latlon_coords(getvar(wrf_ncfile, "RAINC"))
        track_lon.append(lons[np.where(slp == slp.min())].values[0][0])
        track_lat.append(lats[np.where(slp == slp.min())].values[0][0])

    mslp_timeseries = xr.concat(slp_min, dim="time")
    ws_timeseries = xr.concat(ws_max, dim="time")
    return {
        "track_lon": track_lon,
        "track_lat": track_lat,
        "mslp": mslp_timeseries,
        "ws": ws_timeseries,
    }


def plot_track_intensity(harvey, track_details_array, labels, colors):

    fig, axs = plt.subplots(2, 1, figsize=(8, 7))

    axs[0].plot(harvey["date"], harvey["vmax"], "k-", label="OBS")

    for track_details, labs, cols in zip(track_details_array, labels, colors):
        axs[0].plot(
            track_details["ws"]["Time"], track_details["ws"] * 1.95, cols, label=labs
        )

#    axs[0].plot([harvey["date"][42], harvey["date"][42]], [10, 150], "b--")
    axs[0].set_xlabel("Date")
    axs[0].set_ylabel("10 m Wind Speed (knots)")
    axs[0].set_xlim((harvey["date"][33], harvey["date"][50]))
#    axs[0].set_ylim((10, 150))
    myFmt = mdates.DateFormatter("%dZ%H")
    axs[0].xaxis.set_major_formatter(myFmt)

    axs[1].plot(harvey["date"], harvey["mslp"], "k-", label="OBS")
    for track_details, labs, cols in zip(track_details_array, labels, colors):
        axs[1].plot(
            track_details["mslp"]["Time"], track_details["mslp"], cols, label=labs
        )

    plt.legend()
#    axs[1].plot([harvey["date"][42], harvey["date"][42]], [880, 1008], "b--")
    axs[1].set_xlabel("Date")
    axs[1].set_ylabel("MSLP (hPa)")
    axs[1].set_xlim((harvey["date"][33], harvey["date"][50]))
#    axs[1].set_ylim((880, 1008))
    myFmt = mdates.DateFormatter("%dZ%H")
    axs[1].xaxis.set_major_formatter(myFmt)
    plt.tight_layout()


def plot_track(harvey, track_details_array, labels, colors):
    shapefile_path = (
        "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    )
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    coast.plot_coast(ax)
    for geometry in geometries:
        ax.add_geometries(
            [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
        )
    ax.plot(harvey["lon"], harvey["lat"], "k-", label="OBS")

    # 	ax.plot(track_lon[6:], track_lat[6:], 'r-', label='3 km (9-3 km moving)')

    for track_details, labs, cols in zip(track_details_array, labels, colors):
        ax.plot(
            track_details["track_lon"], track_details["track_lat"], cols, label=labs
        )

    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    ax.set_xlim([-102, -83])
    ax.set_ylim([18, 35])
    plt.legend(loc="upper left")

    plt.tight_layout()


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))


home_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_FNL/'

wrf_runs = home_folder + "WRFV4_mp06_cu05/"
wrfoutfile1 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = home_folder + "WRFV4_mp06_cu05/"
wrfoutfile2 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = home_folder + "WRFV4_mp06_cu05/"
wrfoutfile3 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))


out1 = get_track_details(wrfoutfile1)
out2 = get_track_details(wrfoutfile2)
out3 = get_track_details(wrfoutfile3)

plot_track_intensity(harvey, (out1, out2, out3), labels=("MP06", "MP08", 'MP10'), colors=("r", "b", "g"))
plot_track(harvey, (out1, out2, out3), labels=("MP06", "MP08", 'MP10'), colors=("r", "b", 'g'))

plt.show()






