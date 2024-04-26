from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
#from self_utils import ahps, coast
import src.coast as coast
import numpy as np
from scipy.interpolate import griddata
from src.wrf_src import wrf_assign_coords
import pandas as pd
from matplotlib.colors import LogNorm

plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


def rmse(predictions, targets):
    return np.sqrt(np.nanmean((predictions - targets) ** 2))


def correlation_without_nans(array1, array2):
    # Remove NaN values from both arrays
    valid_indices = np.logical_and(~np.isnan(array1), ~np.isnan(array2))
    valid_array1 = array1[valid_indices]
    valid_array2 = array2[valid_indices]

    # Calculate correlation coefficient
    correlation = np.corrcoef(valid_array1, valid_array2)[0, 1]

    return correlation




def get_insitu_pcp(wrfoutfile, labels='LULC 2000'):
 wrf_pcp = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") + getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (getvar(Dataset(wrfoutfile[0]), "RAINC") + getvar(Dataset(wrfoutfile[0]), "RAINNC"))
 slp = getvar(Dataset(wrfoutfile[0]), "slp")

 wrf_lat, wrf_lon = latlon_coords(slp)
 wrf_pcp = wrf_assign_coords(wrf_pcp)

 insitu = pd.read_csv(
    f"../../WPC_pcp_obs/p24i_201708{date}_sortbyvalue.txt", skiprows=3, delimiter="\s+"
)

 domain_bb = [-97.02, -94, 28.11, 30.67]

 insitu = insitu[
    (insitu['Longitude'] >= domain_bb[0]) &
    (insitu['Longitude'] <= domain_bb[1]) &
    (insitu['Latitude'] >= domain_bb[2]) &
    (insitu['Latitude'] <= domain_bb[3])
]

 insitu_pcp_locs = insitu["Precipitation"]
 wrf_pcp_locs = griddata(
    (wrf_lon.values.ravel(), wrf_lat.values.ravel()),
    wrf_pcp.values.ravel(),
    ((insitu["Longitude"], insitu["Latitude"])),
)

# dataset = (insitu_pcp_locs.values*25.4, wrf_pcp_locs)
 dataset = (wrf_pcp_locs - (insitu_pcp_locs.values*25.4) )
 rms = rmse(wrf_pcp_locs, insitu_pcp_locs.values*25.4)
 cor = correlation_without_nans(insitu_pcp_locs.values*25.4, wrf_pcp_locs)
 #title = f"201708{date} | LULC 2000 | RMSE: {np.round(rms, 2)}, R: {np.round(cor, 2)}"
 title = f"{labels} |  {np.round(rms, 2)} / {np.round(cor, 2)}"
 
 return insitu, dataset, title


domain_bb = [-100, -93, 25.5, 31.5]
domain_bb = [-97.52, -94, 28.11, 30.67]
levels = np.array((0.1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500))

num_colors = 8
cmap = plt.get_cmap("bwr", num_colors)

date = '27'

get_wrffiles = lambda base_path: (sorted(glob.glob(base_path + f'/wrfout_d02_2017-*{int(date)}*'))) 


wrfoutfile_2001 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_cntl/test/em_real/')
wrfoutfile_2017 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF_cntl/test/em_real/')
wrfoutfile_2050 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2050/WRF/test/em_real/')
wrfoutfile_2100 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2100/WRF/test/em_real/')

insitu_2001, dataset_2001, title_2001 = get_insitu_pcp(wrfoutfile_2001, labels='(a) LULC 2000') 
insitu_2017, dataset_2017, title_2017 = get_insitu_pcp(wrfoutfile_2017, labels='(b) LULC 2017') 
insitu_2050, dataset_2050, title_2050 = get_insitu_pcp(wrfoutfile_2050, labels='(c) LULC 2050') 
insitu_2100, dataset_2100, title_2100 = get_insitu_pcp(wrfoutfile_2100, labels='(d) LULC 2100') 

fig, ax = plt.subplots(
    2,
    2,
    figsize=(12, 7),
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)
plt.subplots_adjust(left=0.25, right=0.95, wspace=0.2, hspace=0.1)

cont = ax[0, 0].scatter(insitu_2001["Longitude"].values,  insitu_2001["Latitude"].values,  100,  c=dataset_2001,  cmap=cmap, vmin=-400, vmax=400 )
cont = ax[0, 1].scatter(insitu_2017["Longitude"].values,  insitu_2017["Latitude"].values,  100,  c=dataset_2017,  cmap=cmap, vmin=-400, vmax=400 )
cont = ax[1, 0].scatter(insitu_2050["Longitude"].values,  insitu_2050["Latitude"].values,  100,  c=dataset_2050,  cmap=cmap, vmin=-400, vmax=400 )
cont = ax[1, 1].scatter(insitu_2100["Longitude"].values,  insitu_2100["Latitude"].values,  100,  c=dataset_2100,  cmap=cmap, vmin=-400, vmax=400 )

gl = [coast.plot_coast(ax, houston=True) for ax in ax.ravel()]
gl[1].left_labels=False
gl[3].left_labels=False
#    ax.set_xlim((domain_bb[0], domain_bb[1]))
#    ax.set_ylim((domain_bb[2], domain_bb[3]))
label = (title_2001, title_2017, title_2050, title_2100)
[ax.ravel()[index].set_title(label[index]) for index in (0, 1, 2, 3)]

plt.tight_layout()


cbar = plt.colorbar(cont, ax=ax.ravel(), shrink=0.8)
cbar.set_label('Rainfall Error (mm)', fontsize=14) 

plt.savefig(f"../figures_draft01/extended_fig02.jpeg")
#plt.show()







