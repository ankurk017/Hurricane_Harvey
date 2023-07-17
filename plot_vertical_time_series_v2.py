import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
from self_utils import coast
from src.wrf_src import wrf_assign_coords, area_mean
import progressbar
from wrf import to_np, getvar, CoordPair, vertcross, interplevel
import glob
from matplotlib.patches import Rectangle
import numpy as np
import xarray as xr
import matplotlib

# matplotlib.use('Agg')




plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"

wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre_UCM/WRF/test/em_real/wrfout_d02_2017-*")
)[:36]
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF/test/em_real/wrfout_d02_2017-*")
)[:36]

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

pre_wrffiles = wrfoutfile_pre[:index]
post_wrffiles = wrfoutfile_post[:index]


urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.1}

pre_w_profile = []
post_w_profile = []


for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles, post_wrffiles)):
    ncfile_pre = Dataset(prefiles)
    pres = getvar(ncfile_pre, "pres") / 100
    w_pre = wrf_assign_coords(getvar(ncfile_pre, "wa", units="m/s"))

    ncfile_post = Dataset(postfiles)
    pres = getvar(ncfile_post, "pres") / 100
    w_post = wrf_assign_coords(getvar(ncfile_post, "wa", units="m/s"))

    pre_w_profile.append(area_mean(
        interplevel(w_pre, pres, np.arange(1000, 700, -5)), urban_change
    ))
    post_w_profile.append(area_mean(
        interplevel(w_post, pres, np.arange(1000, 700, -5)), urban_change
    ))

fig, axs = plt.subplots(2, 1, figsize=(10, 12))
xr.concat(pre_w_profile, dim='Time').T.plot(ax=axs[0], levels=np.arange(-1.2, 1.4, 0.1))
xr.concat(post_w_profile, dim='Time').T.plot(ax=axs[1], levels=np.arange(-1.2, 1.4, 0.1))
axs[0].set_title('LULC 2001')
axs[1].set_title('LULC 2017')
[ax.invert_yaxis() for ax in axs]
plt.tight_layout()
plt.savefig('../figures/w_wind/w_wind.jpeg')
fig, axs = plt.subplots(1, 1, figsize=(12, 6))
(xr.concat(post_w_profile, dim='Time')-xr.concat(pre_w_profile, dim='Time')).T.plot(ax=axs, levels=np.arange(-0.6, 0.8, 0.1))
axs.invert_yaxis()
plt.tight_layout()
plt.savefig('../figures/w_wind/w_wind_diff1.jpeg')

fig, axs = plt.subplots(1, 1, figsize=(12, 6))
(xr.concat(post_w_profile, dim='Time')-xr.concat(pre_w_profile, dim='Time')).T.plot.contourf(ax=axs, levels=np.arange(-0.6, 0.8, 0.1))
axs.invert_yaxis()
plt.tight_layout()
plt.savefig('../figures/w_wind/w_wind_diff2.jpeg')

#plt.show()


fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
wrf_assign_coords(getvar(Dataset(pre_wrffiles[14]), "wspd_wdir10")).sel(wspd_wdir='wspd').plot(cmap='binary', levels=np.arange(0, 30, 2))

# Define the coordinates of the rectangle
south = urban_change["south_north"] - urban_change["box"] / 2
north = urban_change["south_north"] + urban_change["box"] / 2
west = urban_change["west_east"] - urban_change["box"] / 2
east = urban_change["west_east"] + urban_change["box"] / 2

# Create a rectangle patch
rectangle = Rectangle((west, south), urban_change["box"], urban_change["box"], facecolor='red', alpha=0.5)

# Add the rectangle to the axes
ax.add_patch(rectangle)


coast.plot_coast(ax)
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()
for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')

gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([-96.66, -93.73])
ax.set_ylim([28.21, 30.93])

plt.tight_layout()
plt.savefig('../figures/w_wind/region_gpm_mean_cross.jpeg')
#plt.show()



