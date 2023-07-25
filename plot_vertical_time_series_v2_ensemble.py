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

case = '_cntl'

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"

wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre_UCM/WRF{case}/test/em_real/wrfout_d02_2017-*")
)[:36]
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF{case}/test/em_real/wrfout_d02_2017-*")
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

w_wind_diff = xr.concat(post_w_profile, dim='Time')-xr.concat(pre_w_profile, dim='Time')

fig, axs = plt.subplots(1, 1, figsize=(8.5, 5))
w_wind_diff.T.plot(ax=axs, levels=np.arange(-0.6, 0.8, 0.1), cbar_kwargs={'label': "Vertical Wind Speed (m/s)"})
axs.invert_yaxis()
axs.set_ylabel('Pressure Levels (hPa)')
#axs.xaxis.set_major_locator(mdates.HourLocator(interval=4))

plt.title(f'{case.upper()[1:]}')
plt.tight_layout()
plt.savefig(f'../figures/w_wind/w_wind_diff1{case}.jpeg')

#plt.show()


