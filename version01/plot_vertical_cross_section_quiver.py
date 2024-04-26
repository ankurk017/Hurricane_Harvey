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
from src.wrf_src import wrf_assign_coords
import progressbar
from wrf import to_np, getvar, CoordPair, vertcross
import glob

import matplotlib
matplotlib.use('Agg')

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/"

start_point = CoordPair(lat=29.3, lon=-96.7)
end_point = CoordPair(lat=30.2, lon=-93.9)


pre_wrffiles = sorted(
    glob.glob(
        home_2512
        + f"/WRF_mov1_GFS_IC25_12UTC_v2/pre/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_*"
    )
)
post_wrffiles = sorted(
    glob.glob(
        home_2512
        + f"/WRF_mov1_GFS_IC25_12UTC_v2/post/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_*"
    )
)


home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_2dom/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

pre_wrffiles = wrfoutfile_pre[:index]
post_wrffiles = wrfoutfile_post[:index]




def get_uv(prenc,  postnc):
 pre_u_cross = vertcross(wspd, getvar(prenc, "ua", units="m/s"), wrfin=prenc, start_point=start_point, end_point=end_point, latlon=True, meta=True,)
 pre_v_cross = vertcross(wspd, getvar(prenc, "va", units="m/s"), wrfin=prenc, start_point=start_point, end_point=end_point, latlon=True, meta=True,)
 post_u_cross = vertcross(wspd, getvar(postnc, "ua", units="m/s"), wrfin=postnc, start_point=start_point, end_point=end_point, latlon=True, meta=True,)
 post_v_cross = vertcross(wspd, getvar(postnc, "va", units="m/s"), wrfin=postnc, start_point=start_point, end_point=end_point, latlon=True, meta=True,)
 return pre_u_cross, pre_v_cross, post_u_cross, post_v_cross



for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles, post_wrffiles)):
    ncfile_pre = Dataset(prefiles)
    z = getvar(ncfile_pre, "z")
    wspd = getvar(ncfile_pre, "wa", units="m/s")
    pre_wspd_cross = vertcross(
        wspd,
        z,
        wrfin=ncfile_pre,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )

    ncfile_post = Dataset(postfiles)
    z = getvar(ncfile_post, "z")
    wspd = getvar(ncfile_post, "wa", units="m/s")
    post_wspd_cross = vertcross(
        wspd,
        z,
        wrfin=ncfile_post,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )
    u1, v1, u2, v2 = get_uv(ncfile_pre, ncfile_post)
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    wspd_contours = ax.contourf(
        to_np(post_wspd_cross) - to_np(pre_wspd_cross),
        levels=np.arange(-5, 5.5, 0.5),
        cmap=get_cmap("bwr"),
    )
    ax.quiver(u1.values, v1.values, color='b')
    ax.quiver(u2.values, v2.values, color='r')
    plt.colorbar(wspd_contours, ax=ax)

    coord_pairs = to_np(pre_wspd_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    ax.set_xticks(x_ticks[::5])
    ax.set_xticklabels(x_labels[::5], rotation=45)

    vert_vals = to_np(pre_wspd_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    ax.set_yticks(v_ticks[::20])
    ax.set_yticklabels(np.round(vert_vals[::20] / 1000, 1))
    # ax.set_ylim((0, 12))
    ax.set_xlabel("Latitude, Longitude")
    ax.set_ylabel("Height (m)")

    plt.title("Vertical Cross Section of vertical wind speed (m/s)")
    plt.tight_layout()
    plt.savefig(f"../figures/cross_section_w_uv/w_{str(wspd.Time.values)[:13]}.jpeg")
    plt.close()
    # plt.show()
#plt.show()
    # PLOT TRACK

fig = plt.figure(figsize=(8, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
wrf_assign_coords(getvar(ncfile_post, "wspd_wdir10")).sel(wspd_wdir='wspd').plot(cmap='binary', levels=np.arange(0, 30, 2))
ax.plot((end_point.lon, start_point.lon), (end_point.lat, start_point.lat), 'r-')
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
ax.set_xlim([-102, -92.5])
ax.set_ylim([25.5, 32.7])

plt.tight_layout()
plt.savefig('../figures/cross_section_w_uv/cross.jpeg')

