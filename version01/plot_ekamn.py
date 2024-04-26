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


urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.2}
urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.1}

for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles, post_wrffiles)):
    ncfile_pre = Dataset(prefiles)
    pres = getvar(ncfile_pre, "pres") / 100
    u_pre = wrf_assign_coords(getvar(ncfile_pre, "ua", units="m/s"))
    v_pre = wrf_assign_coords(getvar(ncfile_pre, "va", units="m/s"))
    w_pre = wrf_assign_coords(getvar(ncfile_pre, "wa", units="m/s"))

    ncfile_post = Dataset(postfiles)
    pres = getvar(ncfile_post, "pres") / 100
    u_post = wrf_assign_coords(getvar(ncfile_post, "ua", units="m/s"))
    v_post = wrf_assign_coords(getvar(ncfile_post, "va", units="m/s"))
    w_post = wrf_assign_coords(getvar(ncfile_post, "wa", units="m/s"))

    pre_u_profile = area_mean(
        interplevel(u_pre, pres, np.arange(1000, 700, -5)), urban_change
    )
    pre_v_profile = area_mean(
        interplevel(v_pre, pres, np.arange(1000, 700, -5)), urban_change
    )
    pre_w_profile = area_mean(
        interplevel(w_pre, pres, np.arange(1000, 700, -5)), urban_change
    )
    post_u_profile = area_mean(
        interplevel(u_post, pres, np.arange(1000, 700, -5)), urban_change
    )
    post_v_profile = area_mean(
        interplevel(v_post, pres, np.arange(1000, 700, -5)), urban_change
    )
    post_w_profile = area_mean(
        interplevel(w_post, pres, np.arange(1000, 700, -5)), urban_change
    )

    plt.figure(figsize=(15, 8))

    ax1 = plt.subplot(131)
    ax1_w = plt.subplot(132)
    ax3 = plt.subplot(233)

    ax1.plot(
        np.sqrt(pre_u_profile**2 + pre_v_profile**2),
        pre_u_profile["level"],
        "b",
        label="LULC 2001",
    )
    ax1.plot(
        np.sqrt(post_u_profile**2 + post_v_profile**2),
        post_u_profile["level"],
        "r",
        label="LULC 2017",
    )
    ax1.legend()
    ax1.set_xlabel("Wind Speed (m/s)")
    ax1.set_ylabel("Pressure Levels (hPa)")
    ax1.invert_yaxis()
    ax1.grid(True)
    ax1.set_title("Wind Speed")

    ax1_w.plot(pre_w_profile**2, pre_w_profile["level"], "b", label="LULC 2017")
    ax1_w.plot(post_w_profile**2, post_w_profile["level"], "r", label="LULC 2017")
    ax1_w.set_xlabel("Vertical Wind Speed (m/s)")
    ax1_w.invert_yaxis()
    ax1_w.grid(True)

    ax1_w.set_title("Vertical Wind Speed ")

    ax3.plot(pre_u_profile, pre_v_profile, "b")
    ax3.plot(post_u_profile, post_v_profile, "r")
    ax3.set_xlabel("u wind (m/s)")
    ax3.set_ylabel("v wind (m/s)")
    ax3.grid(True)
    plt.suptitle(str(u_pre.Time.values)[:13], fontweight="bold")

    plt.tight_layout()
    plt.savefig(f"../figures/wind_profiles/profiles_{str(u_pre.Time.values)[:13]}.jpeg")
    plt.close()
#    plt.show()


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
plt.savefig('../figures/wind_profiles/region_vertical_wind_profiles.jpeg')
#plt.show()


