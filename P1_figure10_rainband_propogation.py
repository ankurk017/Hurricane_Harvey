from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from wrf import (
    getvar,
)
import glob
import pandas as pd
import progressbar
from src.wrf_src import find_common_min_max
import itertools
import xarray as xr
from src.wrf_src import wrf_assign_coords, plot_crossline
from wrf import CoordPair

import tropycal.tracks as tracks
from P1_get_rainbands_locs import get_rainbands_locs_updated

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

def get_rainband_cross_section(input_rainfall, ref_var):
 output_rainband = []

 #for rainband in ('rainband1', 'rainband2'):
 for rainband in ('rainband3', 'rainband4'):
    start_point = (get_rainbands_locs_updated()[rainband]['start'].lon, get_rainbands_locs_updated()[rainband]['start'].lat)
    end_point = (get_rainbands_locs_updated()[rainband]['end'].lon, get_rainbands_locs_updated()[rainband]['end'].lat)

    lats = np.linspace(start_point[1], end_point[1], 100)
    lons = np.linspace(start_point[0], end_point[0], 100)

    interpolated_rainfall = [input_rainfall.interp(west_east=lons[index], south_north=lats[index], method='cubic') for index in range(lons.shape[0])]
    cross_sec = xr.concat(interpolated_rainfall, dim='west_east')
    cross_sec['time'] = ref_var["Time"].values
    output_rainband.append(cross_sec)

 return output_rainband




def get_precip(wrf_runs, var_name = 'precip', dates=['27', '28'], var_level=12):

    wrfoutfile = [sorted(glob.glob(wrf_runs + f"wrfout_d02*-{date}_*")) for date in dates]
    wrfoutfile = list(itertools.chain.from_iterable(wrfoutfile))

    rainband_output = []

    for timeid in progressbar.progressbar(range(len(wrfoutfile) - 1)):

        wrf_ncfile2 = Dataset(wrfoutfile[timeid + 1])
        wrf_ncfile1 = Dataset(wrfoutfile[timeid])
        ref_var = getvar(wrf_ncfile1, "uvmet10_wspd_wdir")
        if var_name == 'precip':
            var_rainfall = wrf_assign_coords((
            getvar(wrf_ncfile2, "RAINC") + getvar(wrf_ncfile2, "RAINNC")
        ) - (getvar(wrf_ncfile1, "RAINC") + getvar(wrf_ncfile1, "RAINNC")))
        else:
            var_rainfall = wrf_assign_coords(getvar(wrf_ncfile2, var_name)) 
        var_rainfall = var_rainfall.isel(bottom_top=var_level) if len(var_rainfall.shape)==3 else var_rainfall
        rainband_output.append(get_rainband_cross_section(var_rainfall, ref_var))
    return rainband_output

pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
variable = 'precip'

rainband_pre = get_precip(pre_files, dates=['27', '28', ], var_name = variable)
rainband_post = get_precip(post_files, dates=['27', '28', ], var_name = variable)

pre_rainband1 = xr.concat([rainband_pre[index][0] for index in range(len(rainband_pre))], dim='time')
pre_rainband2 = xr.concat([rainband_pre[index][1] for index in range(len(rainband_pre))], dim='time')

post_rainband1 = xr.concat([rainband_post[index][0] for index in range(len(rainband_pre))], dim='time')
post_rainband2 = xr.concat([rainband_post[index][1] for index in range(len(rainband_pre))], dim='time')

pre_rainband1 = pre_rainband1.where(pre_rainband1 >= 30, 0)
pre_rainband2 = pre_rainband2.where(pre_rainband2 >= 30, 0)
post_rainband1 = post_rainband1.where(post_rainband1 >= 30, 0)
post_rainband2 = post_rainband2.where(post_rainband2 >= 30, 0)


fig, axs = plt.subplots(1, 2, figsize=(13, 7), sharex=True, sharey=True)

minmax = find_common_min_max((pre_rainband1.values, post_rainband1.values))
levels = np.linspace(10, minmax[1]*1.2, 15)
pre_rainband1.plot(ax=axs[0], levels=levels, cmap='Reds')
post_rainband1.plot(ax=axs[1], levels=levels, cmap='Reds')
#plt.savefig('../figures_paper/rainband/pcp_rainband1_loc2_rainband34.jpeg', dpi=300)

fig, axs = plt.subplots(1, 2, figsize=(13, 7), sharex=True, sharey=True)
minmax = find_common_min_max((pre_rainband2.values, post_rainband2.values))
levels = np.linspace(10, minmax[1]*1.2, 15)
pre_rainband2.plot(ax=axs[0], levels=levels, cmap='Reds')
post_rainband2.plot(ax=axs[1], levels=levels, cmap='Reds')
#plt.savefig('../figures_paper/rainband/pcp_rainband2_loc2_rainband34.jpeg', dpi=300)


plot_crossline(glob.glob(f'{pre_files}/wrfout*2017-08-27*')[-5], get_rainbands_locs_updated(), cmap='GnBu', linewidth=4, rainband=('rainband3', 'rainband4'))
#plt.savefig('../figures_paper/rainband/rainband_locs_loc2_rainband34.jpeg', dpi=300)
#plt.show()




variable = 'wa'

rainband_pre = get_precip(pre_files, dates=['27', '28', ], var_name = variable)
rainband_post = get_precip(post_files, dates=['27', '28', ], var_name = variable)

pre_rainband1 = xr.concat([rainband_pre[index][0] for index in range(len(rainband_pre))], dim='time')
pre_rainband2 = xr.concat([rainband_pre[index][1] for index in range(len(rainband_pre))], dim='time')
post_rainband1 = xr.concat([rainband_post[index][0] for index in range(len(rainband_pre))], dim='time')
post_rainband2 = xr.concat([rainband_post[index][1] for index in range(len(rainband_pre))], dim='time')

#pre_rainband1 = pre_rainband1.where(pre_rainband1 >= 30, 0)
#pre_rainband2 = pre_rainband2.where(pre_rainband2 >= 30, 0)
#post_rainband1 = post_rainband1.where(post_rainband1 >= 30, 0)
#post_rainband2 = post_rainband2.where(post_rainband2 >= 30, 0)

fig, axs = plt.subplots(1, 2, figsize=(13, 7), sharex=True, sharey=True)

minmax = find_common_min_max((pre_rainband1.values, post_rainband1.values))
levels = np.linspace(minmax[0], minmax[1]*1.2, 15)
pre_rainband1.plot(ax=axs[0], levels=levels, cmap='bwr')
post_rainband1.plot(ax=axs[1], levels=levels, cmap='bwr')
plt.savefig('../figures_paper/rainband/pcp_rainband1_loc2_rainband34.jpeg', dpi=300)

fig, axs = plt.subplots(1, 2, figsize=(13, 7), sharex=True, sharey=True)
minmax = find_common_min_max((pre_rainband2.values, post_rainband2.values))
levels = np.linspace(minmax[0], minmax[1]*1.2, 15)
pre_rainband2.plot(ax=axs[0], levels=levels, cmap='bwr')
post_rainband2.plot(ax=axs[1], levels=levels, cmap='bwr')
#plt.savefig('../figures_paper/rainband/pcp_rainband2_loc2_rainband34.jpeg', dpi=300)


plot_crossline(glob.glob(f'{pre_files}/wrfout*2017-08-27*')[-5], get_rainbands_locs_updated(), cmap='GnBu', linewidth=4, rainband=('rainband3', 'rainband4'))
#plt.savefig('../figures_paper/rainband/rainband_locs_loc2_rainband34.jpeg', dpi=300)
plt.show()











