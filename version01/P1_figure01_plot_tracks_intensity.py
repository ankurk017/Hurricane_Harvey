import matplotlib.pyplot as plt
import glob
from src.wrf_src import wrf_assign_coords

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
wrf_pbl01 = sorted(glob.glob(f'{home}/wrfout_d01_2017-*'))[::6]

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
wrf_pbl02 = sorted(glob.glob(f'{home}/wrfout_d01_2017-*'))[::6]



out_pbl1 = src.wrf_track.get_track_details(wrf_pbl01)
out_pbl2 = src.wrf_track.get_track_details(wrf_pbl02)


out = (out_pbl1, out_pbl2)

lab = ('LULC 2001', 'LULC 2017', )
col = ("b", 'r',)
titles = 'LULC'


intensity_ax = src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
plt.savefig('../figures_paper/track_pre_post.jpeg', dpi=300)
#intensity_ax.set_xlim(())

#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
plt.savefig('../figures_paper/intensity_pre_post.jpeg', dpi=300)
#track_ax.set_xlim((-105, -85))
#track_ax.set_ylim((20, 45))
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
#src.wrf_track.calculate_error(out, harvey)
#pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures_paper/{titles}_track_stats.csv')

plt.show()





############### get size
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)
from netCDF4 import Dataset



def calculate_radius_of_34_knot_winds(wind_speed_matrix, center_x=None, center_y=None, threshold_knots=34):
    if center_x is None:
        center_x = wind_speed_matrix.shape[1] // 2  # Default to the center of the domain
    if center_y is None:
        center_y = wind_speed_matrix.shape[0] // 2  # Default to the center of the domain
    
    # Create a grid of coordinates
    x, y = np.meshgrid(np.arange(wind_speed_matrix.shape[1]), np.arange(wind_speed_matrix.shape[0]))
    
    # Calculate distance from the hurricane center
    distance = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    # Find the maximum distance where wind speed >= threshold_knots
    radius = np.max(distance[wind_speed_matrix >= threshold_knots])
    
    return radius

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
wrffiles_pre = sorted(glob.glob(f'{home}/wrfout_d02_2017-*'))[::6]

radius_34_knots_pre = []
time_value = []
lats, lons = latlon_coords(getvar(ncfile, "RAINC"))

for index in range(len(wrffiles_pre)):
	print(index)
	ncfile = Dataset(wrffiles_pre[index])
	slp = getvar(ncfile, 'slp')
	track_lon = lons[np.where(slp == slp.min())].values[0][0]
	track_lat = lats[np.where(slp == slp.min())].values[0][0]

	ws = wrf_assign_coords(np.sqrt((getvar(ncfile, "uvmet10", units='knots').sel(u_v='u')**2 + getvar(ncfile, "uvmet10", units='knots').sel(u_v='v')**2 )))

	wind_speed_matrix  = ws.sel(south_north=slice(track_lat-2, track_lat+2), west_east=slice(track_lon-2, track_lon+2))
	radius_34_knots_pre.append(calculate_radius_of_34_knot_winds(wind_speed_matrix)*3)
	time_value.append(slp.Time)

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
wrffiles_post = sorted(glob.glob(f'{home}/wrfout_d02_2017-*'))[::6]

radius_34_knots_post = []
time_value = []
lats, lons = latlon_coords(getvar(ncfile, "RAINC"))

for index in range(len(wrffiles_post)):
        print(index)
        ncfile = Dataset(wrffiles_post[index])
        slp = getvar(ncfile, 'slp')
        track_lon = lons[np.where(slp == slp.min())].values[0][0]
        track_lat = lats[np.where(slp == slp.min())].values[0][0]
        
        ws = wrf_assign_coords(np.sqrt((getvar(ncfile, "uvmet10", units='knots').sel(u_v='u')**2 + getvar(ncfile, "uvmet10", units='knots').sel(u_v='v')**2 )))
        wind_speed_matrix  = ws.sel(south_north=slice(track_lat-2, track_lat+2), west_east=slice(track_lon-2, track_lon+2))

        radius_34_knots_post.append(calculate_radius_of_34_knot_winds(wind_speed_matrix)*3)
        time_value.append(slp.Time)



plt.plot(xr.concat(time_value, dim='Time'), radius_34_knots_pre, 'b', label='LULC 2001')
plt.plot(xr.concat(time_value, dim='Time'), radius_34_knots_post, 'r', label='LULC 2017')
plt.legend()
plt.show()



# plotting contour of 34 kt winds for each timestep
for pre_files, post_files in zip(wrffiles_pre, wrffiles_post):
        ncfile = Dataset(pre_files)
        slp = getvar(ncfile, 'slp')
        track_lon = lons[np.where(slp == slp.min())].values[0][0]
        track_lat = lats[np.where(slp == slp.min())].values[0][0]
        
        ws = wrf_assign_coords(np.sqrt((getvar(ncfile, "uvmet10", units='knots').sel(u_v='u')**2 + getvar(ncfile, "uvmet10", units='knots').sel(u_v='v')**2 )))

        wind_speed_matrix_pre  = ws.sel(south_north=slice(track_lat-2, track_lat+2), west_east=slice(track_lon-2, track_lon+2))

        ncfile = Dataset(post_files)
        slp = getvar(ncfile, 'slp')
        track_lon = lons[np.where(slp == slp.min())].values[0][0]
        track_lat = lats[np.where(slp == slp.min())].values[0][0]
        
        ws = wrf_assign_coords(np.sqrt((getvar(ncfile, "uvmet10", units='knots').sel(u_v='u')**2 + getvar(ncfile, "uvmet10", units='knots').sel(u_v='v')**2 )))
        wind_speed_matrix_post  = ws.sel(south_north=slice(track_lat-2, track_lat+2), west_east=slice(track_lon-2, track_lon+2))

	fig, ax = plt.subplots(1, 2, figsize=(10, 6))
	ax[0].contourf(wind_speed_matrix_pre)
	ax[0].contour(wind_speed_matrix_pre, levels = (32,))
	ax[1].contourf(wind_speed_matrix_post)
	ax[1].contour(wind_speed_matrix_post, levels = (32,))
	plt.close()
plt.show()


