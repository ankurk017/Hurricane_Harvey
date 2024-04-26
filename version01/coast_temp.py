import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def plot_coast(axes: cartopy.mpl.geoaxes.GeoAxes, color='black', linewidth=2,  gridlines_alpha=0.5) -> None:
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

    return gl


