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
import cartopy.feature as cfeature


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})



basin = tracks.TrackDataset(basin='north_atlantic',
                            source='ibtracs', include_btk=False)
harvey = basin.get_storm(('harvey', 2017))


# PLOT TRACK
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()

file_path = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06/test.txt'

A=pd.read_csv(file_path, header=None, delimiter=r"\s+")


fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
#coast.plot_coast(ax)

ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)


for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')
ax.plot(harvey['lon'], harvey['lat'], 'k-', label='OBS')
ax.plot(A[3],A[2],'r-',transform=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([-102, -83])  
ax.set_ylim([18, 35]) 
plt.legend(loc='upper left')

plt.tight_layout()
#plt.savefig('../figures/May22_Harvey_track.jpeg')

from datetime import datetime
df = pd.read_csv(file_path, delimiter='\s+', header=None, names=['Type', 'DateTime', 'Lat', 'Lon', 'Pressure', 'WS'])

df['DateTime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d_%H:%M:%S')
df.set_index('DateTime', inplace=True)

harvey_dt = harvey.to_dataframe()
harvey_dt.set_index('date', inplace=True)

plt.figure(figsize=(11, 7))
df['Pressure'].plot(color='r')
harvey_dt['mslp'].plot(color='k')
plt.ylabel('Minimum Sea Level Pressure (hPa)')

plt.figure(figsize=(11, 7))
df['WS'].plot(color='r')
harvey_dt['vmax'].plot(color='k')
plt.ylabel('Minimum Sea Level Pressure (hPa)')
plt.show()





