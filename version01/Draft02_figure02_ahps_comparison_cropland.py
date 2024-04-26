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
from matplotlib import gridspec
from src.wrf_src import wrf_assign_coords
import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd

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



def plot_ahps(date, plot_axis):
	ahps_home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"
	ahps_files = glob.glob(ahps_home + "nws_precip*_conus.nc")
	ahps_fileid = np.where([f'201708{date}' in files for files in ahps_files])[0][0]
	ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])

	cont = plot_axis.contourf(ahps_lon, ahps_lat, ahps_pcp,
			      levels=levels, 
			      locator=ticker.LogLocator(), 
			      cmap='gist_ncar',) #extend='both')

	coast.plot_coast(plot_axis, houston=True)
	plot_axis.set_xlim((domain_bb[0], domain_bb[1]))
	plot_axis.set_ylim((domain_bb[2], domain_bb[3]))
	plot_axis.set_title(f'(a) AHPS (201708{date})', fontsize=13.5)
	plot_axis.text(np.nanmax(ahps_lon), np.nanmin(ahps_lat), f'{np.nanmean(ahps_pcp).round(2)} mm', 
               verticalalignment='bottom', horizontalalignment='right', fontsize=14, color='blue')

	return {'ahps_lon': ahps_lon, 'ahps_lat':ahps_lat, 'ahps_pcp':ahps_pcp}


def plot_accp(wrfoutfile, ahps_grids, plot_axis, year_label='2001'):
	ahps_lon = ahps_grids['ahps_lon']
	ahps_lat = ahps_grids['ahps_lat']
	ahps_pcp = ahps_grids['ahps_pcp']
	wrf_pcp = (
	    getvar(Dataset(wrfoutfile[-1]), "RAINC") +
	    getvar(Dataset(wrfoutfile[-1]), "RAINNC")
	) - (
	    getvar(Dataset(wrfoutfile[0]), "RAINC") +
	    getvar(Dataset(wrfoutfile[0]), "RAINNC")
	)
	slp = getvar(Dataset(wrfoutfile[0]), "slp")
	wrf_lat, wrf_lon = latlon_coords(slp)



	cont = plot_axis.contourf(wrf_lon, wrf_lat, wrf_pcp, cmap="gist_ncar",
			      levels=levels, 
			      locator=ticker.LogLocator(), #extend='both'
			      )

	coast.plot_coast(plot_axis, houston=True)


	plot_axis.set_xlim((domain_bb[0], domain_bb[1]))
	plot_axis.set_ylim((domain_bb[2], domain_bb[3]))


	merged_lats = np.arange(domain_bb[2], domain_bb[3], 0.05)
	merged_lons = np.arange(domain_bb[0], domain_bb[1], 0.05)
	merged_meshlon, merged_meshlat = np.meshgrid(merged_lons, merged_lats)


	ahps_pcp_domain = griddata((ahps_lon.ravel(), ahps_lat.ravel(
	)), ahps_pcp.ravel(), (merged_meshlon, merged_meshlat))
	wrf_pcp_domain = griddata((wrf_lon.values.ravel(), wrf_lat.values.ravel(
	)), wrf_pcp.values.ravel(), (merged_meshlon, merged_meshlat))

	rms = rmse(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
	cor = correlation_without_nans(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
	#title = f"{ahps_files[ahps_fileid].split('/')[-1].split('_')[-2]} | {case} | RMSE: {np.round(rms, 2)}, R: {np.round(cor, 2)}"
	plot_axis.set_title(f'{year_label} | {np.round(rms, 2)} / {np.round(cor, 2)}', fontsize=13.25)

	plot_axis.text(np.nanmax(wrf_lon)-1.5, np.nanmin(wrf_lat), f'{np.nanmean(wrf_pcp_domain).round(2)} mm', 
               verticalalignment='bottom', horizontalalignment='right', fontsize=14, color='blue')


	return cont

get_wrffiles = lambda base_path: (sorted(glob.glob(base_path + f'/wrfout_d02_2017-*{int(date)-1}*')) + sorted(glob.glob(base_path + f'/wrfout_d02_2017-*{int(date)}*')))[12:12+25]
#get_wrffiles = lambda base_path: (sorted(glob.glob(base_path + f'/wrfout_d02_2017-*{int(date)}*'))) 
get_allwrffiles = lambda base_path: (sorted(glob.glob(base_path + f'/wrfout_d02_*'))) 

domain_bb = [-100, -93, 25.5, 31.5]
levels = np.array((10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 900))

date = '27'

wrfoutfile_2001 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_cntl/test/em_real/')
wrfoutfile_2017 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF_cntl/test/em_real/')
wrfoutfile_2050 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2050/WRF/test/em_real/')
wrfoutfile_2100 = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2100/WRF/test/em_real/')

wrfoutfile_2050_crop = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2050_cropland/WRF/test/em_real/')
wrfoutfile_2100_crop = get_wrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2100_cropland/WRF/test/em_real/')


fig = plt.figure(figsize=(14.5, 8))
gs = gridspec.GridSpec(2, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.4, left=0.07, right=0.9)

ahps_grid_ax = plt.subplot(gs[0, 0], projection=ccrs.PlateCarree())
cont_2001_ax = plt.subplot(gs[1, 0], projection=ccrs.PlateCarree())
#cont_2017_ax = plt.subplot(gs[1, 0], projection=ccrs.PlateCarree())
cont_2050_ax = plt.subplot(gs[0, 1], projection=ccrs.PlateCarree())
cont_2100_ax = plt.subplot(gs[0, 2], projection=ccrs.PlateCarree())

cont_2050crop_ax = plt.subplot(gs[1, 1], projection=ccrs.PlateCarree())
cont_2100crop_ax = plt.subplot(gs[1, 2], projection=ccrs.PlateCarree())
colorbar_ax = plt.subplot(gs[:, -1])

ahps_grid = plot_ahps(date=date, plot_axis=ahps_grid_ax)
plot_accp(wrfoutfile_2001, ahps_grid, plot_axis=cont_2001_ax, year_label='(d) LULC 2001')
#plot_accp(wrfoutfile_2017, ahps_grid, plot_axis=cont_2017_ax, year_label='(d) LULC 2017')

plot_accp(wrfoutfile_2050, ahps_grid, plot_axis=cont_2050_ax, year_label='(b) LULC 2050')
cont = plot_accp(wrfoutfile_2100, ahps_grid, plot_axis=cont_2100_ax, year_label='(c) LULC 2100')


plot_accp(wrfoutfile_2050_crop, ahps_grid, plot_axis=cont_2050crop_ax, year_label='(e) LULC 2050 Crop')
plot_accp(wrfoutfile_2100_crop, ahps_grid, plot_axis=cont_2100crop_ax, year_label='(f) LULC 2100 Crop')


cbar = plt.colorbar(cont, cax=colorbar_ax, ticks=levels)
cbar.ax.set_yticklabels(levels)
cbar.ax.set_ylabel('Precipitation (mm)')

plt.tight_layout()
plt.savefig('../figures_draft01/fig02a_tmp.jpeg', dpi=400)


# basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
# harvey = basin.get_storm(("harvey", 2017))

# wrfoutfile_2001 = get_allwrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_cntl/test/em_real/')[::3]
# wrfoutfile_2050 = get_allwrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2050/WRF/test/em_real/')[::3]
# wrfoutfile_2100 = get_allwrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2100/WRF/test/em_real/')[::3]
# wrfoutfile_2050_crop = get_allwrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2050_cropland/WRF/test/em_real/')[::3]
# wrfoutfile_2100_crop = get_allwrffiles('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_simulations//WRF_FNL_2100_cropland/WRF/test/em_real/')[::3]


# out_pbl1 = src.wrf_track.get_track_details(wrfoutfile_2001)
# out_pbl3 = src.wrf_track.get_track_details(wrfoutfile_2050)
# out_pbl4 = src.wrf_track.get_track_details(wrfoutfile_2100)
# out_pbl3_cropland = src.wrf_track.get_track_details(wrfoutfile_2050_crop)
# out_pbl4_cropland = src.wrf_track.get_track_details(wrfoutfile_2100_crop)

# out = (out_pbl1, out_pbl3, out_pbl4, out_pbl3_cropland, out_pbl4_cropland)

# lab = ['2001', '2050', '2100', '2050_crop', '2100_crop']
# col = ("b", 'r','g', 'm', '')
# titles = 'EPA_LULC'

# #track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col, extent=[-103, -80, 16, 32])
# track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col, extent=[-98.54, -91.72, 26, 32], figsize=(4.5, 4), legend_fontsize=9)
# track_ax.set_title('(a) Tracks')
# #track_ax.set_xlim((-105, -85))
# #track_ax.set_ylim((20, 45))
# #plt.savefig(f'../figures/{titles}_plot_track.jpeg')

# src.wrf_track.calculate_error(out, harvey)
# pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures_draft01/{titles}_track_stats.csv')

# plt.savefig('../figures_draft01/fig02b_tmp.jpeg', dpi=400)

# #plt.show()


#intensity_ax = src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig('../figures_draft01/intenisty.jpeg', dpi=400)