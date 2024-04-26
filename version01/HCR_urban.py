import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np
from scipy.interpolate import griddata
import cartopy.feature as cfeature

from matplotlib import ticker, cm
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import matplotlib.pyplot as plt
from wrf import getvar, interplevel, to_np, latlon_coords
from src.wrf_src import wrf_assign_coords
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"

case = 'pre'
name = 'LULC 2001'

wrfoutfile = sorted(glob.glob(
    f'/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/{case}/WRF_9-3-51/WRFV4/wrfout_d02_2017-*27*'))


wrf_wa = getvar(Dataset(wrfoutfile[-1]), "wa")
z = getvar(Dataset(wrfoutfile[-1]), "z", units="m")
p = getvar(Dataset(wrfoutfile[-1]), "pressure")
wspd_500 = interplevel(wrf_wa, z, 500)
slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)


bb = np.array((wrf_lon.min().values, wrf_lon.max().values,
              wrf_lat.min().values, wrf_lat.max().values))
bb = [-100, -94.5, 26, 31]
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()


ax = plt.subplot(projection=ccrs.PlateCarree())
wrf_assign_coords(wspd_500).plot(ax=ax, cmap='bwr')

ax.coastlines()
ax.add_feature(cfeature.STATES.with_scale('50m'))

for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                      facecolor='none', edgecolor='blue')

gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=True,
    linewidth=2,
    color="gray",
    alpha=0.5,
    linestyle="--",
)
gl.top_labels = False
gl.bottom_labels = True
gl.right_labels = False
gl.left_labels = True

gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
plt.title(f'w at 350 m AGL (m/s) - {name}')
plt.savefig(f'../figures/HCR_{name}.jpeg')
plt.show()
# plt.close()
