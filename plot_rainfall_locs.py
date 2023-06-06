from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np
from scipy.interpolate import griddata
from src.wrf_src import wrf_assign_coords
import pandas as pd
from matplotlib.colors import LogNorm

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


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


home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"

ahps_files = glob.glob(home + "nws_precip*_conus.nc")
date = 30
case = 'post'

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'
wrfoutfile = sorted(glob.glob(home_2512 + f'/WRF_mov1_GFS_IC25_12UTC_v2/{case}/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_2017-*{date}*'))

wrf_pcp = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") + getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (getvar(Dataset(wrfoutfile[0]), "RAINC") + getvar(Dataset(wrfoutfile[0]), "RAINNC"))
slp = getvar(Dataset(wrfoutfile[0]), "slp")

wrf_lat, wrf_lon = latlon_coords(slp)
wrf_pcp = wrf_assign_coords(wrf_pcp)


insitu = pd.read_csv(
    "../../WPC_pcp_obs/p24i_20170828_sortbyvalue.txt", skiprows=3, delimiter="\s+"
)
insitu_pcp_locs = insitu["Precipitation"]
wrf_pcp_locs = griddata(
    (wrf_lon.values.ravel(), wrf_lat.values.ravel()),
    wrf_pcp.values.ravel(),
    ((insitu["Longitude"], insitu["Latitude"])),
)

dataset = (insitu_pcp_locs.values*25.4, wrf_pcp_locs)
dataset = (insitu_pcp_locs.values*25.4 - wrf_pcp_locs, )

domain_bb = [-100, -93, 25.5, 31.5]
domain_bb = [-97.52, -94, 28.11, 30.67]
levels = np.array((0.1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500))

fig, ax = plt.subplots(
    1,
    1,
    figsize=(9, 5),
    sharey=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)

for id in range(1):
    cont = ax.scatter(
        insitu["Longitude"].values,
        insitu["Latitude"].values,
        80,
        c=dataset,
        cmap="bwr",
#        norm=LogNorm(),
    )

    rms = rmse(wrf_pcp_locs, insitu_pcp_locs.values*25.4)
    cor = correlation_without_nans(insitu_pcp_locs.values*25.4, wrf_pcp_locs)
    title = f"WRF 201708{date} | RMSE: {np.round(rms, 2)}, R: {np.round(cor, 2)}"


    shapefile_path = (
        "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
    )
    reader = shpreader.Reader(shapefile_path)
    geometries = reader.geometries()

    for geometry in geometries:
        ax.add_geometries(
            [geometry], ccrs.PlateCarree(), facecolor="none", edgecolor="blue"
        )

    coast.plot_coast(ax)
    ax.set_xlim((domain_bb[0], domain_bb[1]))
    ax.set_ylim((domain_bb[2], domain_bb[3]))
    ax.set_title(title)
plt.tight_layout()
plt.colorbar(cont)
plt.savefig(f"../figures/{case}_WRF_pcp_locs_201708{date}.jpeg")
#plt.show()

