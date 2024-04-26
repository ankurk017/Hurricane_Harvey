from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

from wrf import (
    to_np,
    getvar,
    smooth2d,
    get_cartopy,
    cartopy_xlim,
    cartopy_ylim,
    latlon_coords,
    interplevel,
)

# from self_utils import ahps, coast
import src.coast as coast
import src.ahps as ahps
import numpy as np
from scipy.interpolate import griddata
from matplotlib import gridspec
from src.wrf_src import wrf_assign_coords
import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
import xarray as xr
import metpy
from metpy.interpolate import cross_section
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import progressbar
from src.wrf_src import wrf_assign_coords, plot_bb, area_mean

plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "font.family": "monospace"}
)


wrf_files_pre = sorted(
    glob.glob(
        f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_HR/test/em_real/wrfout_d01_2017-08-27*"
    )
) #[6:13]
wrf_files_post = sorted(
    glob.glob(
        f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_HR/test/em_real/wrfout_d01_2017-08-27*"
    )
) #[6:13]

wrf_files_pre = sorted(glob.glob('/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*'))[:70]
wrf_files_post = sorted(glob.glob('/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*'))[:70]

urban_change = {"south_north": 29.74, "west_east": -95.49, "box":0.7}

outdir = '../figures_draft01/radial_tangential_1000-800_HR_profile_FNL/'


def process_hr(_wrf_files_post, pres_levels = np.arange(1000, 10, -10)):
    p, hr = [], []
    for index in progressbar.progressbar(range(len(_wrf_files_post))):
        prefix = _wrf_files_post[index].split('/')[-1][11:24]
        ncfile = Dataset(_wrf_files_post[index])
        p.append(getvar(ncfile, "pressure"))
        hr.append(getvar(ncfile, "H_DIABATIC"))
    p_xr = wrf_assign_coords(xr.concat(p, dim='time'))
    hr_xr = wrf_assign_coords(xr.concat(hr, dim='time'))

    p_bb = area_mean(xr.concat(p_xr, dim='time'), urban_change)
    hr_bb = area_mean(xr.concat(hr_xr, dim='time'), urban_change)
    hr_output = xr.concat([hr_bb.isel(time=time_id).assign_coords(bottom_top=p_bb.isel(time=time_id)).interp(bottom_top=pres_levels) for time_id in range(hr_bb.time.shape[0])], dim='time')
    hr_output['time'] = hr_output['Time']
    hr_output = hr_output.rename({'bottom_top':'Pressure'})
    return hr_output


post_hr = process_hr(wrf_files_post, pres_levels=np.linspace(1000, 800, 25))
pre_hr = process_hr(wrf_files_pre, pres_levels=np.linspace(1000, 800, 25))


fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharex=True, sharey=True)
pre_hr.T.plot.contourf(ax=axs[0], levels=np.linspace(0, 0.0020, 20), cmap='jet')
post_hr.T.plot.contourf(ax=axs[1], levels=np.linspace(0, 0.0020, 20), cmap='jet')
axs[0].invert_yaxis()
plt.show()


'''
 fig, axs = plt.subplots(1, 1, figsize=(6, 7))
 axs.plot(hr_pre[-1], hr_pre[-1].level, 'b', marker='o', label='LULC 2001')
 axs.plot(hr_post[-1], hr_post[-1].level, 'r', marker='o', label='LULC 2017')
 axs.set_xlabel('Heating Rate (K/s)')
 axs.invert_yaxis()
 axs.grid(True)
 axs.set_ylabel('Pressure (hPa)')
 plt.legend()
 plt.title(prefix)
 #plt.show()
 plt.savefig(f'{outdir}/{prefix}.jpeg', dpi=400, bbox_inches='tight')
'''



