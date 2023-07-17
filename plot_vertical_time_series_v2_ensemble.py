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



plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

case='_min3'

home_2512 = sorted(glob.glob(f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/*_UCM/WRF{case}*"))
#home_2512 = sorted(glob.glob(f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/*_UCM/WRF"))

fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

for home_folder in home_2512:

 wrfoutfile_pre = sorted(
    glob.glob(home_folder + f"/test/em_real/wrfout_d02_2017-*")
)[:36][12:]


 urban_change = {"south_north": 29.79, "west_east": -95.72, "box":0.1}

 pre_w_profile = []

 for prefiles in progressbar.progressbar(wrfoutfile_pre):
    ncfile_pre = Dataset(prefiles)
    pres = getvar(ncfile_pre, "pres") / 100
    w_pre = wrf_assign_coords(getvar(ncfile_pre, "wa", units="m/s"))

    pre_w_profile.append(area_mean(
        interplevel(w_pre, pres, np.array((950, 850))), urban_change
    ))


 xr.concat(pre_w_profile, dim='Time').sel(level=950).plot(ax=axs[0], label=home_folder.split('/')[-2].split('_')[0])

 xr.concat(pre_w_profile, dim='Time').sel(level=850).plot(ax=axs[1], label=home_folder.split('/')[-2].split('_')[0])
plt.legend()
plt.suptitle(case)
plt.tight_layout()
plt.savefig(f'../figures/w_wind/ensemble{case}.jpeg')
plt.show()


