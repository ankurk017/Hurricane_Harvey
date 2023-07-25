import numpy as np
from matplotlib.patches import Rectangle
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from wrf import getvar
from netCDF4 import Dataset
from self_utils import ahps, coast


def wrf_assign_coords(var):
    lon_1d = var.XLONG.mean(dim="south_north").values
    lat_1d = var.XLAT.mean(dim="west_east").values
    return var.drop(["XLONG", "XLAT"]).assign_coords(
        south_north=lat_1d, west_east=lon_1d
    )


def renamelatlon(var):
    return var.rename({"west_east": "longitudes", "south_north": "latitudes"})


def area_mean(pre_wspd, location, box=0.2):
    lat_id = np.where(
        np.logical_and(
            pre_wspd["south_north"].values > location["south_north"] - location["box"],
            pre_wspd["south_north"].values <= location["south_north"] + location["box"],
        )
    )[0]
    lon_id = np.where(
        np.logical_and(
            pre_wspd["west_east"].values > location["west_east"] - location["box"],
            pre_wspd["west_east"].values <= location["west_east"] + location["box"],
        )
    )[0]

    pre_wspd_mean = (
        pre_wspd.isel(south_north=lat_id, west_east=lon_id)
        .mean(dim="south_north")
        .mean(dim="west_east")
    )
    return pre_wspd_mean


def crop_region(pre_wspd, location, box=0.2):
    lat_id = np.where(
        np.logical_and(
            pre_wspd["south_north"].values > location["south_north"] - location["box"],
            pre_wspd["south_north"].values <= location["south_north"] + location["box"],
        )
    )[0]
    lon_id = np.where(
        np.logical_and(
            pre_wspd["west_east"].values > location["west_east"] - location["box"],
            pre_wspd["west_east"].values <= location["west_east"] + location["box"],
        )
    )[0]

    crop_region = pre_wspd.isel(south_north=lat_id, west_east=lon_id)
    return crop_region


def plot_bb(wrf_file, location):
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    wrf_assign_coords(getvar(Dataset(wrf_file), "wspd_wdir10")).sel(
        wspd_wdir="wspd"
    ).plot(cmap="binary", levels=np.arange(0, 30, 2))

    south = location["south_north"] - location["box"] / 2
    north = location["south_north"] + location["box"] / 2
    west = location["west_east"] - location["box"] / 2
    east = location["west_east"] + location["box"] / 2

    rectangle = Rectangle(
        (west, south),
        location["box"],
        location["box"],
        facecolor="red",
        alpha=0.5,
    )

    ax.add_patch(rectangle)

    coast.plot_coast(ax)
    shapefile_path = (
        "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    )
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()
    for geometry in geometries:
        ax.add_geometries(
            [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
        )

    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    ax.set_xlim([-96.66, -93.73])
    ax.set_ylim([28.21, 30.93])

    plt.tight_layout()

    return fig


def plot_crossline(wrf_file, start_point, end_point):
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    wrf_assign_coords(getvar(Dataset(wrf_file), "wspd_wdir10")).sel(wspd_wdir="wspd").plot(
        cmap="binary", levels=np.arange(0, 30, 2)
    )
    ax.plot((end_point.lon, start_point.lon), (end_point.lat, start_point.lat), "r-")
    coast.plot_coast(ax)
    shapefile_path = (
        "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    )
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()
    for geometry in geometries:
        ax.add_geometries(
            [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
        )

    gl = ax.gridlines(draw_labels=True)
    gl.right_labels = False
    gl.top_labels = False
    ax.set_xlim([-97.5, -94.2])
    ax.set_ylim([29, 31])

    plt.tight_layout()
    return fig










