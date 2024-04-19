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


case = '_cntl'

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
pre_wrffiles = sorted(
    glob.glob(home_2512 + f"/pre/WRF{case}/test/em_real/wrfout_d02_2017-08-26*")
)
post_wrffiles = sorted(
    glob.glob(home_2512 + f"/post/WRF{case}/test/em_real/wrfout_d02_2017-08-26*")
)



location = (-95.499,  29.74)  # Houston location
box = 0.25

crop_lonlat_id = {'west_east':slice(location[0]-box, location[0]+box), 'south_north':slice(location[1]-box, location[1]+box)}
heights = np.arange(10, 1000, 10)
heights = (0, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000)

UU, VV = [], []
for index in progressbar.progressbar(range(len(post_wrffiles))):
 ncfile_pre = Dataset(pre_wrffiles[index])
 u = wrf_assign_coords(getvar(ncfile_pre, "ua")).sel(crop_lonlat_id)
 v = wrf_assign_coords(getvar(ncfile_pre, "va")).sel(crop_lonlat_id)

 z = wrf_assign_coords(getvar(ncfile_pre, "z", units="m")).sel(crop_lonlat_id)
 
 UU.append(interplevel(u, z, heights).mean(dim=['south_north', 'west_east']))
 VV.append(interplevel(v, z, heights).mean(dim=['south_north', 'west_east']))

U = xr.concat(UU, dim='Time').mean(dim='Time')
V = xr.concat(VV, dim='Time').mean(dim='Time')


fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

axes[0].plot(U, V)
[axes[0].text(U[index], V[index], str(heights[index])) for index in range(len(heights))]
axes[0].axhline(0, color='blue', linestyle='--', linewidth=1.5)
axes[0].axvline(0, color='blue', linestyle='--', linewidth=1.5)
axes[0].grid()
axes[0].set_xlabel('u (m/s)')
axes[0].set_ylabel('v (m/s)')

axes[1].quiver(np.ones(U.shape[0]), heights, U, V, scale=20, width=.02, color='blue')
axes[1].set_xlabel('Wind Vectors (m/s)')
axes[1].set_ylabel('Height (m)')
axes[1].set_title('Vertical Profile of Wind Vectors')
axes[1].grid()

#axes[1].set_aspect(0.0001, adjustable='box')

plt.tight_layout()

plt.show()


















