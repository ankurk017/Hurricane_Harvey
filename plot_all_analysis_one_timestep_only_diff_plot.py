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
from src.wrf_src import wrf_assign_coords, plot_crossline
import progressbar
from wrf import to_np, getvar, CoordPair, vertcross
import glob
from plot_rainfall_func import plot_rainfall_and_winds
import matplotlib

#matplotlib.use("Agg")

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def get_uv(prenc, postnc):
    pre_u_cross = vertcross(
        getvar(prenc, "ua", units="m/s"),
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=prenc,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )
    pre_v_cross = vertcross(
        getvar(prenc, "wa", units="m/s"),
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=prenc,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )
    post_u_cross = vertcross(
        getvar(postnc, "ua", units="m/s"),
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=postnc,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )
    post_v_cross = vertcross(
        getvar(postnc, "wa", units="m/s"),
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=postnc,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )
    return pre_u_cross, pre_v_cross * 10, post_u_cross, post_v_cross * 10


start_point = CoordPair(lat=29.3, lon=-96.5)
end_point = CoordPair(lat=30.5, lon=-94.5)

start_point = CoordPair(lat=29.55, lon=-96.07)
end_point = CoordPair(lat=30.27, lon=-94.84)

case = '_cntl'

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"
pre_wrffiles = sorted(
    glob.glob(home_2512 + f"/pre_UCM/WRF{case}/test/em_real/wrfout_d02_2017-*")
)[21:24]# [22:24]
post_wrffiles = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF{case}/test/em_real/wrfout_d02_2017-*")
)[21:24] #[22:24]

for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles[1:], post_wrffiles[1:])):
    ncfile_pre = Dataset(prefiles)
    pres = getvar(ncfile_pre, "pres") / 100
    wspd = getvar(ncfile_pre, "wa", units="m/s")
    wspd = getvar(ncfile_pre, "QVAPOR")
    pre_wspd_cross = vertcross(
        wspd,
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=ncfile_pre,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )

    ncfile_post = Dataset(postfiles)
    pres = getvar(ncfile_post, "pres") / 100
    wspd = getvar(ncfile_post, "wa", units="m/s")
    wspd = getvar(ncfile_post, "QVAPOR")
    post_wspd_cross = vertcross(
        wspd,
        pres,
        levels=np.arange(1000, 700, -5),
        wrfin=ncfile_post,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    )

    u1, v1, u2, v2 = get_uv(ncfile_pre, ncfile_post)

    coord_pairs = to_np(pre_wspd_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    vert_vals = to_np(pre_wspd_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])

    fig, ax = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    cbar=plt.colorbar(wspd_contours, ax=ax[1], orientation='horizontal', pad=0.25)
    cbar.set_label('Water Vapor Mixing Ratio (kg/kg)')

    wspd_contours = ax[2].contourf(
        x_ticks,
        vert_vals,
        post_wspd_cross,
        #levels=np.arange(-6, 6.5, 0.5), cmap=get_cmap("bwr")
        #levels=np.arange(0, 170, 20),
        cmap=get_cmap("coolwarm"),
        extend="both",
    )
    ax[2].quiver(x_ticks[::4], vert_vals[::4], u2[::4, ::4]-u1[::4, ::4], v2[::4, ::4]-v1[::4, ::4])
    cbar=plt.colorbar(wspd_contours, ax=ax[2], orientation='horizontal', pad=0.25)
    cbar.set_label('Water Vapor Mixing Ratio Difference (kg/kg)')
    for axs in ax:
        
        axs.set_xticks(x_ticks[::5])
        axs.set_xticklabels(x_labels[::5], rotation=45)
        axs.invert_yaxis()
        axs.set_xlabel("Latitude, Longitude")
        axs.set_ylabel("Pressure (hPa)")

    ax[0].set_title(f"LULC 2001: {str(wspd['Time'].values)[:13]}")
    ax[1].set_title(f"LULC 2017: {str(wspd['Time'].values)[:13]}")
    plt.tight_layout()

    #    plt.savefig(f"../figures/cross_section_cwlp_uw_actual_UCM_ensemble/{case}_w_{str(wspd.Time.values)[:13]}.jpeg")
    #    plt.close()

# plt.show()
# PLOT TRACK

plt.show()






