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
import src.ahps as ahps
import numpy as np
from scipy.interpolate import griddata


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

date = '28'
case = 'WRF_PBL04'
#case = 'post'

labels = {"WRF_PBL00":'No PBL', "WRF_PBL01":'YSU', "WRF_PBL02":'MYJ TKE', "WRF_PBL04":'EDMF QNSE'}


yyear = '2001' if case=='pre' else '2017'

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"

ahps_files = glob.glob(home + "nws_precip*_conus.nc")
home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/PBL/'
wrfoutfile = sorted(glob.glob(home_2512 + f'/{case}/test/em_real/wrfout_d02_2017-*{date}*'))
wrfoutfile = (sorted(glob.glob(home_2512 + f'/{case}/test/em_real/wrfout_d02_2017-*{int(date)-1}*')) + sorted(glob.glob(home_2512 + f'/{case}/test/em_real/wrfout_d02_2017-*{int(date)}*')))[12:12+25]

wrf_pcp = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") +
    getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (
    getvar(Dataset(wrfoutfile[0]), "RAINC") +
    getvar(Dataset(wrfoutfile[0]), "RAINNC")
)
slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)

ahps_fileid = np.where([f'201708{date}' in files for files in ahps_files])[0][0]
ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])


domain_bb = [-100, -93, 25.5, 31.5]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900))

fig, ax = plt.subplots(
    1, 2, figsize=(14, 5), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)
cont = ax[0].contourf(ahps_lon, ahps_lat, ahps_pcp,
                      levels=levels, 
                      locator=ticker.LogLocator(), 
                      cmap='gist_ncar',) # extend='both')

shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()

for geometry in geometries:
    ax[0].add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')


coast.plot_coast(ax[0])
ax[0].set_xlim((domain_bb[0], domain_bb[1]))
ax[0].set_ylim((domain_bb[2], domain_bb[3]))
ax[0].set_title('AHPS')

cont = ax[1].contourf(wrf_lon, wrf_lat, wrf_pcp, cmap="gist_ncar",
                      levels=levels, 
                      locator=ticker.LogLocator(), #extend='both'
                      )

shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()

for geometry in geometries:
    ax[1].add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')

coast.plot_coast(ax[1])


ax[1].set_xlim((domain_bb[0], domain_bb[1]))
ax[1].set_ylim((domain_bb[2], domain_bb[3]))
ax[1].set_title(f'WRF ({labels[case]})')
cbar = plt.colorbar(cont, ax=ax.ravel(), ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')


merged_lats = np.arange(domain_bb[2], domain_bb[3], 0.05)
merged_lons = np.arange(domain_bb[0], domain_bb[1], 0.05)
merged_meshlon, merged_meshlat = np.meshgrid(merged_lons, merged_lats)


ahps_pcp_domain = griddata((ahps_lon.ravel(), ahps_lat.ravel(
)), ahps_pcp.ravel(), (merged_meshlon, merged_meshlat))
wrf_pcp_domain = griddata((wrf_lon.values.ravel(), wrf_lat.values.ravel(
)), wrf_pcp.values.ravel(), (merged_meshlon, merged_meshlat))

rms = rmse(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
cor = correlation_without_nans(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
#print(rms, cor)
title = f"{ahps_files[ahps_fileid].split('/')[-1].split('_')[-2]} | {labels[case]} | RMSE: {np.round(rms, 2)}, R: {np.round(cor, 2)}"

plt.savefig(f"../figures_paper/ahps/12UTC_2km_AHPS_{case}_pcp_{ahps_files[ahps_fileid].split('/')[-1].split('_')[-2]}.jpeg")
#plt.show()

########## plotting differenbce
levels = np.arange(-500, 500, 100)
fig, ax = plt.subplots(
    1, 1, figsize=(8, 5), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
)
cont = ax.contourf(merged_lons, merged_lats, wrf_pcp_domain-ahps_pcp_domain,
                      levels=levels, cmap='bwr', extend='both')

shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()

for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')

coast.plot_coast(ax)


ax.set_xlim((domain_bb[0], domain_bb[1]))
ax.set_ylim((domain_bb[2], domain_bb[3]))
cbar = plt.colorbar(cont, ax=ax, ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation Error (mm)')
title = f"{ahps_files[ahps_fileid].split('/')[-1].split('_')[-2]} | {labels[case]} | RMSE: {np.round(rms, 2)}, R: {np.round(cor, 2)}"
plt.title(title)
plt.savefig(f"../figures_paper/ahps/12UTC_Diff_2km_AHPS_{case}_pcp_{ahps_files[ahps_fileid].split('/')[-1].split('_')[-2]}.jpeg")
#plt.show()





