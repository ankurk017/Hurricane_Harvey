from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from wrf import (
    getvar, interplevel
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
    wrfoutfile = list(itertools.chain.from_iterable(wrfoutfile))[6:13]

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

        if len(var_rainfall.shape)==3:
            print(f'{var_name} is a 3D variable')
            pres = getvar(wrf_ncfile1, "pres")
            var_rainfall = interplevel(var_rainfall, pres, 750)
        else:
            var_rainfall =  var_rainfall

        rainband_output.append(get_rainband_cross_section(var_rainfall, ref_var))
    return rainband_output

pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"


variables = ['helicity', 'updraft_helicity', 'pvo', 'pw']
labels = ['', '', '', '']
fig, ax = plt.subplots(2, 2, figsize=(13, 8), sharex=True)

for i, var_name in enumerate(variables):
    rainband_pre = get_precip(pre_files, dates=['27'], var_name=var_name)
    rainband_post = get_precip(post_files, dates=['27'], var_name=var_name)
    pre_rainband = xr.concat([rainband_pre[index][0] for index in range(len(rainband_pre))], dim='time')
    post_rainband = xr.concat([rainband_post[index][0] for index in range(len(rainband_pre))], dim='time')

    row, col = divmod(i, 2)  # Calculate row and column for subplot
    pre_rainband.sel(time=slice('2017-08-27T06', '2017-08-27T12')).mean(dim='time').plot(color='b', label='LULC 2001', ax=ax[row, col])
    post_rainband.sel(time=slice('2017-08-27T06', '2017-08-27T12')).mean(dim='time').plot(color='r', label='LULC 2017', ax=ax[row, col])
#    title = f"{post_rainband.attrs['description']} ({post_rainband.attrs['units']})"
#    ax[row, col].set_title(f'{title}')
#    ax[row, col].set_ylabel(labels[i])
    
ax[row, col].legend()

plt.tight_layout()  # To prevent overlapping subplots
plt.show()





