import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
from matplotlib import colors


def get_lulc_colormap():
    lulc_classes = {
        "Evergreen Needleleaf Forest": "#05450a",
        "Evergreen Broadleaf Forest": "#086a10",
        "Deciduous Needleleaf Forest": "#54a708",
        "Deciduous Broadleaf Forest": "#78d203",
        "Mixed Forest": "#009900",
        "Closed Shrublands": "#c6b044",
        "Open Shrublands": "#dcd159",
        "Woody Savannas": "#dade48",
        "Savannas": "#fbff13",
        "Grasslands": "#b6ff05",
        "Permanent wetlands": "#27ff87",
        "Croplands": "#006400",
        "Urban and Built-Up": "#FF0000",
        "cropland/natural vegetation mosaic": "#ADFF2F",
        "Snow and Ice": "#69fff8",
        "Barren or Sparsely Vegetated": "#f9ffa4",
        "Water": "#1c0dff",
    }

    lulc_colormap = ListedColormap(list(lulc_classes.values()))
    return lulc_colormap, lulc_classes
