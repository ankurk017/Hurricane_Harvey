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
import src.wrf_track

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))


home_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_FNL/'

wrf_runs = home_folder + "WRFV4_mp10_cu01/"
wrfoutfile1 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = home_folder + "WRFV4_mp10_cu05/"
wrfoutfile2 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = home_folder + "WRFV4_mp10_cu06/"
wrfoutfile3 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))


out1 = src.wrf_track.get_track_details(wrfoutfile1)
out2 = src.wrf_track.get_track_details(wrfoutfile2)
out3 = src.wrf_track.get_track_details(wrfoutfile3)

src.wrf_track.plot_track_intensity(harvey, (out1, out2, out3), labels=("MP06", "MP08", 'MP10'), colors=("r", "b", "g"))
src.wrf_track.plot_track(harvey, (out1, out2, out3), labels=("MP06", "MP08", 'MP10'), colors=("r", "b", 'g'))

plt.show()






