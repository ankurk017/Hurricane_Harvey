import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


# def plot_coast(axes):
def plot_coast(axes: cartopy.mpl.geoaxes.GeoAxes, color='black', linewidth=2, houston=False, houston_color='blue', houston_linewidth=1.5, gridlines_alpha=0.5) -> None:
    """
    Plot natural features and gridlines on a map using Cartopy.

    Parameters
    ----------
    axes : cartopy.mpl.geoaxes.GeoAxes
        The axes object to plot on.

    Returns
    -------
    None

    """
    countries = cfeature.NaturalEarthFeature(
        scale="10m", category="cultural", name="admin_0_countries", facecolor="none"
    )
    states = cfeature.NaturalEarthFeature(
        scale="10m",
        category="cultural",
        name="admin_1_states_provinces_lines",
        facecolor="none",
    )
    axes.add_feature(countries, edgecolor=color, linewidth=linewidth)
    axes.add_feature(states, edgecolor=color, linewidth=linewidth)
    if houston:
        plot_houston(axes, color=houston_color, linewidth=houston_linewidth)
    gl = axes.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=2,
        color="gray",
        alpha=gridlines_alpha,
        linestyle="--",
    )
    gl.top_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.right_labels = False
    gl.xlines = True
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER




import cartopy.io.shapereader as shpreader

def plot_houston(ax, color='black', linewidth=1):
    shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()
    for geometry in geometries:
             ax.add_geometries([geometry], ccrs.PlateCarree(),
                                 facecolor='none', edgecolor=color, linewidth=linewidth)
    return None

