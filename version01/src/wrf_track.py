#from self_utils import coast
import src.coast  as coast
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
import pandas as pd
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

    fig, axs = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    axs[0].plot(harvey["time"], harvey["vmax"], "k-", label="OBS", marker='o')

    for track_details, labs, cols in zip(track_details_array, labels, colors):
        axs[0].plot(
            track_details["ws"]["Time"], track_details["ws"] * 1.95, cols, label=labs, marker='o'
        )

    #    axs[0].plot([harvey["date"][42], harvey["date"][42]], [10, 150], "b--")
    axs[0].set_xlabel("Date")
    axs[0].set_ylabel("10 m Wind Speed (knots)")
    axs[0].set_xlim((harvey["time"][33], harvey["time"][-15]))
    #    axs[0].set_ylim((10, 150))
    myFmt = mdates.DateFormatter("%dZ%H")
    axs[0].xaxis.set_major_formatter(myFmt)

    axs[1].plot(harvey["time"], harvey["mslp"], color="k", label="OBS", marker='o')
    for track_details, labs, cols in zip(track_details_array, labels, colors):
        axs[1].plot(
            track_details["mslp"]["Time"], track_details["mslp"], color=cols, label=labs, marker='o'
        )

    plt.legend(fontsize=13)
    #    axs[1].plot([harvey["date"][42], harvey["date"][42]], [880, 1008], "b--")
    axs[1].set_xlabel("Date")
    axs[1].set_ylabel("MSLP (hPa)")
    axs[1].set_xlim((harvey["time"][33], harvey["time"][-15]))
    #    axs[1].set_ylim((880, 1008))
    myFmt = mdates.DateFormatter("%dZ%H")
    axs[1].xaxis.set_major_formatter(myFmt)

    [ax.tick_params(axis='x', rotation=30)  for ax in axs]
    [ax.grid(True) for ax in axs]

    plt.tight_layout()
    return axs

def plot_track(harvey, track_details_array, labels, colors, extent=None, figsize=(9,7), legend_fontsize=14):
    shapefile_path = (
        "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    )
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()

    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(left=0.29)

    ax = plt.axes(projection=ccrs.PlateCarree())
    coast.plot_coast(ax)
    for geometry in geometries:
        ax.add_geometries(
            [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
        )
    ax.plot(harvey["lon"], harvey["lat"], "k-", label="OBS", marker='o')

    hurdat_crop = harvey

    for index in np.arange(0, harvey.time.shape[0], 4):
        ax.text(hurdat_crop["lon"][index]+0.05, hurdat_crop["lat"][index], str(hurdat_crop["time"][index])[:13].replace(' ', 'T')[8:], transform=ccrs.PlateCarree(), fontsize=legend_fontsize, color='#7D7D7D')

    for track_details, labs, cols in zip(track_details_array, labels, colors):
        ax.plot(
            track_details["track_lon"], track_details["track_lat"], cols, label=labs, 
        )

    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False

    #    ax.set_xlim([-102, -87.5])
    #    ax.set_ylim([22.5, 35])
    if extent is None:
        ax.set_xlim([-99, -92])
        ax.set_ylim([26, 31])
    else:
        ax.set_xlim([extent[0], extent[1]])
        ax.set_ylim([extent[2], extent[3]])
        
    plt.legend(loc="upper left", fontsize=legend_fontsize)

    plt.tight_layout()
    return ax

def calculate_error(out, harvey):
    ws_error = []
    mslp_error = []
    track_error = []
    for out1a in out:
        model_ws = np.array(
            [
                harvey["vmax"][
                    np.where(
                        harvey["time"]
                        == pd.to_datetime(str(val.values)).to_pydatetime()
                    )[0]
                ]
                for val in out1a["mslp"]["Time"]
            ],
            dtype="object",
        )
        ws_error.append(
            (
                model_ws[np.array([ws.shape[0] != 0 for ws in model_ws])]
                - (out1a["ws"][np.array([ws.shape[0] != 0 for ws in model_ws])]).values
                * 1.95
            ).mean()
        )

        model_mslp = np.array(
            [
                harvey["mslp"][
                    np.where(
                        harvey["time"]
                        == pd.to_datetime(str(val.values)).to_pydatetime()
                    )[0]
                ]
                for val in out1a["mslp"]["Time"]
            ],
            dtype="object",
        )
        mslp_error.append(
            (
                model_mslp[np.array([ws.shape[0] != 0 for ws in model_ws])]
                - (
                    out1a["mslp"][np.array([ws.shape[0] != 0 for ws in model_ws])]
                ).values
            ).mean()
        )

        model_lon = np.array(
            [
                harvey["lon"][
                    np.where(
                        harvey["time"]
                        == pd.to_datetime(str(val.values)).to_pydatetime()
                    )[0]
                ]
                for val in out1a["mslp"]["Time"]
            ],
            dtype="object",
        )
        model2_lon = model_lon[np.array([ws.shape[0] != 0 for ws in model_ws])] - (
            np.array(out1a["track_lon"])[
                np.array([ws.shape[0] != 0 for ws in model_ws])
            ]
        )
        model_lat = np.array(
            [
                harvey["lat"][
                    np.where(
                        harvey["time"]
                        == pd.to_datetime(str(val.values)).to_pydatetime()
                    )[0]
                ]
                for val in out1a["mslp"]["Time"]
            ],
            dtype="object",
        )
        model2_lat = model_lat[np.array([ws.shape[0] != 0 for ws in model_ws])] - (
            np.array(out1a["track_lat"])[
                np.array([ws.shape[0] != 0 for ws in model_ws])
            ]
        )

        track_error.append(
            np.sqrt((np.concatenate(model2_lon**2 + model2_lat**2)).astype('float')).mean() * 111.11
        )

    return {"ws": ws_error, "mslp": mslp_error, "track": track_error}
