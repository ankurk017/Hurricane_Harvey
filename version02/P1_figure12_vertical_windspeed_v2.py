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
from src.coast import plot_coast
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


home_2512 = "/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/"

wrfoutfile_pre = sorted(glob.glob(home_2512 + f"/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-2*00"))[36:-23]
wrfoutfile_post = sorted(glob.glob(home_2512 + f"/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-2*00"))[36:-23]


pre_wrffiles = wrfoutfile_pre
post_wrffiles = wrfoutfile_post

urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.1} # perfectly working

urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.2} # perfectly working
urban_change = {"south_north": 29.75, "west_east": -95.35, "box":0.5} # no urban change

case = 'no-Urban'

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
        interplevel(w_pre, pres, np.arange(1000, 800, -5)), urban_change
    ))
    post_w_profile.append(area_mean(
        interplevel(w_post, pres, np.arange(1000, 800, -5)), urban_change
    ))

fig, axs = plt.subplots(2, 1, figsize=(8, 9), sharex=True)

levels = np.arange(-0.8, 0.9, 0.1)
pre_cont = xr.concat(pre_w_profile, dim='Time').T.plot(ax=axs[0], levels=levels, add_colorbar=False)
xr.concat(post_w_profile, dim='Time').T.plot(ax=axs[1], levels=levels, add_colorbar=False)

axs[0].set_title('LULC 2001')
axs[1].set_title('LULC 2017')

cbar = plt.colorbar(pre_cont, ax=axs.ravel())
cbar.set_label('Vertical Wind (m/s)')
[ax.invert_yaxis() for ax in axs]
[ax.set_ylabel('Pressure (hPa)') for ax in axs]

plt.savefig(f'../../figures_draft02/wwind/{case}_wwind_domain_plot.jpeg', dpi=400, bbox_inches='tight')
#plt.show()



fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
wrf_assign_coords(getvar(Dataset(pre_wrffiles[14]), "RAINNC")).plot(cmap='Blues', levels=np.arange(0, 1000, 50))

# Define the coordinates of the rectangle
south = urban_change["south_north"] - urban_change["box"] / 2
north = urban_change["south_north"] + urban_change["box"] / 2
west = urban_change["west_east"] - urban_change["box"] / 2
east = urban_change["west_east"] + urban_change["box"] / 2

# Create a rectangle patch
rectangle = Rectangle((west, south), urban_change["box"], urban_change["box"], facecolor='red', alpha=0.5)

# Add the rectangle to the axes
ax.add_patch(rectangle)


plot_coast(ax, houston=True)
gl = ax.gridlines(draw_labels=True)
gl.right_labels = False
gl.top_labels = False
ax.set_xlim([-96.66, -93.73])
ax.set_ylim([28.21, 30.93])

plt.tight_layout()
plt.savefig(f'../../figures_draft02/wwind/{case}_wwind_domain.jpeg')
#plt.show()



