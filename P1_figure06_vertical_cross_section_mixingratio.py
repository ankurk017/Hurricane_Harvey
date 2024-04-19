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
from src.wrf_src import wrf_assign_coords, plot_crossline
import progressbar
from wrf import to_np, getvar, CoordPair, vertcross
import glob
from plot_rainfall_func import plot_rainfall_and_winds
import matplotlib
import imageio
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

start_point = CoordPair(lat=29.25, lon=-96.14) # for 2512
end_point = CoordPair(lat=30.27, lon=-94.7)  # for 2512
from P1_get_rainbands_locs import get_rainbands_locs_updated

rainband = 'rainband3'
start_point1 = CoordPair(lat=29.9, lon=-95.87) # for 2512
end_point1 = CoordPair(lat=29.9, lon=-94.22)  # for 2512


case = '_cntl'

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
pre_wrffiles = sorted(
    glob.glob(home_2512 + f"/pre/WRF{case}/test/em_real/wrfout_d03_2017-08-27*")
)[7:15]# [22:24]
post_wrffiles = sorted(
    glob.glob(home_2512 + f"/post/WRF{case}/test/em_real/wrfout_d03_2017-08-27*")
)[7:15] #[22:24]
output_dir = f"../figures_paper/Ensemble/Cross-Section/{case[1:]}/"
image_files = []
#ax = plot_rainfall_and_winds(pre=pre_wrffiles, post=post_wrffiles, plot_cross=True, start_point=start_point, end_point=end_point)
#plt.show()

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

    wspd_contours = ax[0].contourf(
        x_ticks,
        vert_vals,
        pre_wspd_cross,
        #levels=np.arange(-6, 6.5, 0.5), cmap=get_cmap("bwr")
        levels=np.arange(0.01, 0.025, 0.002),
        cmap=get_cmap("rainbow"),
        extend="both",
    )
    ax[0].quiver(x_ticks[::4], vert_vals[::4], u1[::4, ::4], v1[::4, ::4])
    cbar=plt.colorbar(wspd_contours, ax=ax[0], orientation='horizontal', pad=0.25)
    cbar.set_label('Water Vapor Mixing Ratio (kg/kg)')
    wspd_contours = ax[1].contourf(
        x_ticks,
        vert_vals,
        post_wspd_cross,
        #levels=np.arange(-6, 6.5, 0.5), cmap=get_cmap("bwr")
        levels=np.arange(0.01, 0.025, 0.002),
        cmap=get_cmap("rainbow"),
        extend="both",
    )
    ax[1].quiver(x_ticks[::4], vert_vals[::4], u2[::4, ::4], v2[::4, ::4])
    cbar=plt.colorbar(wspd_contours, ax=ax[1], orientation='horizontal', pad=0.25)
    cbar.set_label('Water Vapor Mixing Ratio (kg/kg)')
    cbar.ax.tick_params(rotation=15)

    wspd_contours = ax[2].contourf(
        x_ticks,
        vert_vals,
        post_wspd_cross-pre_wspd_cross,
        #levels=np.arange(-6, 6.5, 0.5), cmap=get_cmap("bwr")
        #levels=np.arange(0, 170, 20),
        cmap=get_cmap("coolwarm"),
        extend="both",
    )
    ax[2].quiver(x_ticks[::4], vert_vals[::4], u2[::4, ::4]-u1[::4, ::4], v2[::4, ::4]-v1[::4, ::4])
    cbar=plt.colorbar(wspd_contours, ax=ax[2], orientation='horizontal', pad=0.25)
    cbar.set_label('Water Vapor Mixing Ratio Difference (kg/kg)')
    cbar.ax.tick_params(rotation=15)
    for axs in ax:
        
        axs.set_xticks(x_ticks[::12])
        axs.set_xticklabels(x_labels[::12], rotation=45, ha='right')
        axs.invert_yaxis()
    #    axs.set_xlabel("Latitude, Longitude")
        axs.set_ylabel("Pressure (hPa)")

    ax[0].set_title(f"LULC 2001 | {str(wspd['Time'].values)[:13]}")
    ax[1].set_title(f"LULC 2017 | {str(wspd['Time'].values)[:13]}")
    ax[2].set_title(f"{case[1:]} | LULC 2017 - LULC 2001: {str(wspd['Time'].values)[:13]}")
    plt.tight_layout()

    image_file = f"{output_dir}/{case}_w_{str(wspd.Time.values)[:13]}.jpeg"
    plt.savefig(image_file)
    image_files.append(image_file)

    plt.close()

# plt.show()
# PLOT TRACK

#plt.show()


## Create a GIF from the saved images
gif_file = output_dir+f"/Cross_Section.gif"
with imageio.get_writer(gif_file, mode="I", duration=0.5) as writer:
    for image_file in image_files:
        image = imageio.imread(image_file)
        writer.append_data(image)


print("GIF created successfully!")

plt.close()
plot_crossline(prefiles, start_point, end_point)
plt.show()




