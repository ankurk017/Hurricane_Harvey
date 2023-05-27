from self_utils import coast
from netCDF4 import Dataset
import matplotlib.dates as mdates
import numpy as np
import matplotlib.pyplot as plt
from wrf import (
    getvar,
    latlon_coords,
)
import cartopy.crs as ccrs
import xarray as xr
import progressbar

import cartopy.io.shapereader as shpreader

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
