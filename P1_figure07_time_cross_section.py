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
import src.coast as coast
from src.wrf_src import wrf_assign_coords, plot_crossline
import progressbar
from wrf import to_np, getvar, CoordPair, vertcross
import glob
from plot_rainfall_func import plot_rainfall_and_winds
import matplotlib
import imageio
#matplotlib.use("Agg")
from P1_get_rainbands_locs import get_rainbands_locs
import xarray as xr

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


rainband = get_rainbands_locs()

rainband1 = rainband['rainband1']
rainband2 = rainband['rainband2']


start_point = rainband2['start']
end_point = rainband2['end']

case = '_cntl'

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
pre_wrffiles = sorted(
    glob.glob(home_2512 + f"/pre/WRF{case}/test/em_real/wrfout_d03_2017-*")
)[36:86]# [22:24]
post_wrffiles = sorted(
    glob.glob(home_2512 + f"/post/WRF{case}/test/em_real/wrfout_d03_2017-*")
)[36:86] #[22:24]
output_dir = f"../figures/Ensemble/Cross-Section/{case[1:]}/"
image_files = []
plot_crossline(pre_wrffiles[-5], rainband1['start'], rainband1['end'])
plot_crossline(pre_wrffiles[-5], rainband2['start'], rainband2['end'])
#plt.show()

pre_wspd_cross, post_wspd_cross =[], []

for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles[1:], post_wrffiles[1:])):
    ncfile_pre = Dataset(prefiles)
    pres = getvar(ncfile_pre, "pres") / 100
    wspd = getvar(ncfile_pre, "wa", units="m/s")
    wspd = getvar(ncfile_pre, "rh")
    pre_wspd_cross.append(vertcross(
        wspd,
        pres,
        levels=(950, 925),
        wrfin=ncfile_pre,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    ).mean(dim='vertical'))

    ncfile_post = Dataset(postfiles)
    pres = getvar(ncfile_post, "pres") / 100
    wspd = getvar(ncfile_post, "wa", units="m/s")
    wspd = getvar(ncfile_post, "rh")
    post_wspd_cross.append(vertcross(
        wspd,
        pres,
        levels=(950, 925),
        wrfin=ncfile_post,
        start_point=start_point,
        end_point=end_point,
        latlon=True,
        meta=True,
    ).mean(dim='vertical'))


fig, axs = plt.subplots(1, 2, figsize=(13, 8), sharey=True)

xr.concat(pre_wspd_cross, dim='Time').plot(ax=axs[0]) #, levels = np.arange(60, 105, 5), extend='both')
xr.concat(post_wspd_cross, dim='Time').plot(ax=axs[1]) # , levels = np.arange(60, 100, 5), extend='both')

plt.tight_layout()
plt.show()



latitudes = np.linspace(start_point.lat, end_point.lat, 20)
longitudes = np.linspace(start_point.lon, end_point.lon, 20)

dbz_cross_pre, dbz_cross_post = [], []

for prefiles, postfiles in progressbar.progressbar(zip(pre_wrffiles[1:], post_wrffiles[1:])):
    ncfile_pre = Dataset(prefiles)
    ncfile_post = Dataset(postfiles)

    dbz_pre = wrf_assign_coords(getvar(ncfile_pre, "mdbz"))
    dbz_post = wrf_assign_coords(getvar(ncfile_post, "mdbz"))

    dbz_cross_pre.append(np.array([dbz_pre.interp(south_north=latitudes[index], west_east=longitudes[index]) for index in range(longitudes.shape[0])]))
    dbz_cross_post.append(np.array([dbz_post.interp(south_north=latitudes[index], west_east=longitudes[index]) for index in range(longitudes.shape[0])]))














