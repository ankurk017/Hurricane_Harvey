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
import progressbar
import matplotlib
matplotlib.use('Agg')

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


domain_bb = [-100, -93, 25.5, 31.5]
domain_bb = [-96.24, -94.57, 29.22, 30.46]

levels = np.array((0.1, 10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700))
levels = np.array((0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200))
levels = np.array((0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2, 3, 4, 5, 7, 10, 15, 20, 40, 50))

#home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'

#prefiles = sorted(glob.glob(home_2512 + f'/WRF_mov1_GFS_IC25_12UTC_v2/pre/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_2017-*'))
#postfiles = sorted(glob.glob(home_2512 + f'/WRF_mov1_GFS_IC25_12UTC_v2/post/WRF_mp10_cu05_no_ocean_physics/wrfout_d01_2017-*'))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_2dom/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

prefiles = wrfoutfile_pre[:index]
postfiles = wrfoutfile_post[:index]

years = (2001, 2017, )
for fileid in progressbar.progressbar(range(len(prefiles)-1)):
	wrf_pcp_pre = getvar(Dataset(prefiles[fileid+1]), 'wspd_wdir10').sel(wspd_wdir='wspd')
	wrf_pcp_post = getvar(Dataset(postfiles[fileid+1]), 'wspd_wdir10').sel(wspd_wdir='wspd')

	slp = getvar(Dataset(prefiles[fileid+1]), "slp")
	wrf_lat, wrf_lon = latlon_coords(slp)

	pre_post = (wrf_pcp_pre, wrf_pcp_post)
	fig, ax = plt.subplots(
	    1, 3, figsize=(17, 4), sharey=True, subplot_kw={"projection": ccrs.PlateCarree()},
	)
	for pre_post_id in range(2):
		cont = ax[pre_post_id].contourf(wrf_lon, wrf_lat, pre_post[pre_post_id], cmap="jet",
				      levels=levels, locator=ticker.LogLocator())
#				      levels=levels, )
		shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
		reader = shpreader.Reader(shapefile_path)
		geometries = reader.geometries()

		for geometry in geometries:
		    ax[pre_post_id].add_geometries([geometry], ccrs.PlateCarree(),
					 facecolor='none', edgecolor='blue')

		coast.plot_coast(ax[pre_post_id])


		ax[pre_post_id].set_xlim((domain_bb[0], domain_bb[1]))
		ax[pre_post_id].set_ylim((domain_bb[2], domain_bb[3]))
		ax[pre_post_id].set_title(f'WRF (LULC {years[pre_post_id]})')
	ax[1].set_yticklabels([])
	cbar = plt.colorbar(cont, ax=ax[:2].ravel(), ticks=levels)
	cbar.ax.set_yticklabels(levels)
#	cbar.ax.set_ylabel('Precipitation (mm)')

	diff_levels = np.arange(-10, 12, 2)
	cont = ax[2].contourf(wrf_lon, wrf_lat, wrf_pcp_post-wrf_pcp_pre, cmap="bwr", levels=diff_levels)
	shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
	reader = shpreader.Reader(shapefile_path)
	geometries = reader.geometries()
	for geometry in geometries:
	    ax[2].add_geometries([geometry], ccrs.PlateCarree(),
				 facecolor='none', edgecolor='blue')
	coast.plot_coast(ax[2])

	ax[2].set_xlim((domain_bb[0], domain_bb[1]))
	ax[2].set_ylim((domain_bb[2], domain_bb[3]))
	ax[2].set_title(f'WRF (LULC {years[pre_post_id]})')
	cbar = plt.colorbar(cont, ax=ax[2], ticks=diff_levels)
	cbar.ax.set_yticklabels(diff_levels)
	cbar.ax.set_ylabel('Wind Speed error (m/s)')

	plt.savefig(f"../figures/houston_ws/Houston_{str(slp.Time.values)[:13]}.jpeg")
	plt.close()

#plt.show()












