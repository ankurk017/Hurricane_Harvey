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
from scipy.interpolate import RegularGridInterpolator
from src.wrf_src import find_common_min_max
import itertools
import xarray as xr

import tropycal.tracks as tracks
from src.wrf_src import wrf_assign_coords


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

def convert_polar1(var_cropped, margin=5, resolution=0.09, coords=['south_north', 'west_east']):
    r = np.arange(0, margin, resolution)
    ang = np.arange(0, 361, 1) * np.pi / 180
    r_mesh, ang_mesh = np.meshgrid(r, ang)

    x_polar = r_mesh * np.cos(ang_mesh)
    y_polar = r_mesh * np.sin(ang_mesh)

    xlong_values = var_cropped[coords[0]]
    xlat_values = var_cropped[coords[1]]

    pcp_values = var_cropped.values

    pcp_interp = RegularGridInterpolator((xlong_values, xlat_values), pcp_values, method='linear', bounds_error=False, fill_value=np.nan)

    x_polar_flat = x_polar.flatten()
    y_polar_flat = y_polar.flatten()

    pcp_polar_flat = pcp_interp((x_polar_flat, y_polar_flat))
    pcp_polar = pcp_polar_flat.reshape(x_polar.shape)


    return r, ang, pcp_polar

def get_precip(wrf_runs, var_name = 'precip', dates=['27', '28'], hurr_center_ref = True, ref_loc = (-97.2, 29.2), obs=None):

    wrfoutfile = [sorted(glob.glob(wrf_runs + f"wrfout_d02*-{date}_*")) for date in dates]
    wrfoutfile = list(itertools.chain.from_iterable(wrfoutfile))
    wrfoutfile = wrfoutfile[13:] if '26' in dates else wrfoutfile
    harvey = obs
    harvey_finer = harvey.to_xarray().interp(
        time=pd.date_range(
            harvey.to_xarray()["time"][0].values,
            harvey.to_xarray()["time"][-1].values,
            freq="1H",
        )
    )

    margin = 4
    wrf_time = []
    precip_polar = []

    for timeid in progressbar.progressbar(range(len(wrfoutfile) - 1)):

        wrf_ncfile2 = Dataset(wrfoutfile[timeid + 1])
        wrf_ncfile1 = Dataset(wrfoutfile[timeid])

        ref_var = getvar(wrf_ncfile1, "uvmet10_wspd_wdir")
        wrf_time.append(ref_var["Time"].values)

        var_rainfall = wrf_assign_coords((
            getvar(wrf_ncfile2, "RAINC") + getvar(wrf_ncfile2, "RAINNC")
        ) - (getvar(wrf_ncfile1, "RAINC") + getvar(wrf_ncfile1, "RAINNC")))

        if hurr_center_ref:
            cen_lon = harvey_finer.sel(time=ref_var["Time"].values)["lon"].values
            cen_lat = harvey_finer.sel(time=ref_var["Time"].values)["lat"].values
        else:
            cen_lon = ref_loc[0]
            cen_lat = ref_loc[1]

        var_cropped = var_rainfall.sel(west_east=slice(cen_lon - margin, cen_lon + margin), 
                               south_north=slice(cen_lat - margin, cen_lat + margin))

        var_cropped = var_cropped.assign_coords(
            west_east=var_cropped.west_east - cen_lon)
        var_cropped = var_cropped.assign_coords(
            south_north=var_cropped.south_north - cen_lat)

        r, ang, var_polar = convert_polar1(var_cropped)
        precip_polar.append(var_polar)
    precip_ds = xr.Dataset(
    data_vars={var_name: (('time', 'angle', 'radius'), np.array(precip_polar))},
    coords={'radius': r*111.11, 'angle': ang, 'time': np.array(wrf_time)})
    return precip_ds


def get_precip_gpm(wrf_runs, var_name = 'precip', dates=['27', '28'], hurr_center_ref = True, ref_loc = (-97.2, 29.2), obs=None):

    harvey = obs
    harvey_finer = harvey.to_xarray().interp(
        time=pd.date_range(
            harvey.to_xarray()["time"][0].values,
            harvey.to_xarray()["time"][-1].values,
            freq="30min",
        )
    )

    margin = 4
    wrf_time = []
    precip_polar = []

    for timeid in progressbar.progressbar(range(wrf_runs.time.shape[0] - 1)):
        var_rainfall = wrf_runs.isel(time=timeid)
        wrf_time.append(var_rainfall['time'].values)

        if hurr_center_ref:
            cen_lon = harvey_finer.sel(time=var_rainfall['time'].values)["lon"].values
            cen_lat = harvey_finer.sel(time=var_rainfall['time'].values)["lat"].values
        else:
            cen_lon = ref_loc[0]
            cen_lat = ref_loc[1]

        var_cropped = var_rainfall.sel(lon=slice(cen_lon - margin, cen_lon + margin),
                               lat=slice(cen_lat - margin, cen_lat + margin))

        var_cropped = var_cropped.assign_coords(
            lon=var_cropped.lon - cen_lon)
        var_cropped = var_cropped.assign_coords(
            lat=var_cropped.lat - cen_lat)

        r, ang, var_polar = convert_polar1(var_cropped, coords=['lon','lat'])
        precip_polar.append(var_polar)
    precip_ds = xr.Dataset(
    data_vars={var_name: (('time', 'angle', 'radius'), np.array(precip_polar))},
    coords={'radius': r*111.11, 'angle': ang, 'time': np.array(wrf_time)})
    return precip_ds




basin = tracks.TrackDataset(basin="north_atlantic", source="hurdat", include_btk=False  )
harvey = basin.get_storm(("harvey", 2017))


pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
variable = 'precip'

pre_precip = get_precip(pre_files, dates=['27', '28', ], var_name = variable, obs=harvey)[variable]
post_precip = get_precip(post_files, dates=['27', '28', ], var_name = variable, obs=harvey)[variable]


pre_precip_polar_mean = pre_precip.mean(dim='time').isel(radius=np.arange(0, 30))
post_precip_polar_mean =  post_precip.mean(dim='time').isel(radius=np.arange(0, 30))

minmax = find_common_min_max((pre_precip_polar_mean.values, post_precip_polar_mean.values))
levels = np.linspace(minmax[0], minmax[1], 15)

pre_precip = get_precip(pre_files, dates=['27', '28'], var_name = variable, hurr_center_ref=True, obs=harvey)[variable]
post_precip = get_precip(post_files, dates=['27', '28'], var_name = variable, hurr_center_ref=True, obs=harvey)[variable]

plot_pre = pre_precip.mean(dim='angle')
plot_post = post_precip.mean(dim='angle')

import datetime


gpm_files = sorted(
    glob.glob("/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/GPM/*nc4")
)[432:750]

gpm = xr.open_mfdataset(gpm_files)["precipitationCal"]
gpm_timestep = [
    datetime.datetime(*dates.timetuple()[:6]) for dates in gpm["time"].values
]
gpm['time'] = gpm_timestep

start_date = '2017-08-27'
end_date = '2017-08-28'

levels = np.arange(0, 16, 1)
cropped_gpm = gpm.sel(time=slice(start_date, end_date))

gpm_polar = get_precip_gpm(cropped_gpm, var_name = 'precip', dates=['27', '28'], hurr_center_ref = True, obs=harvey)['precip'].mean(dim='time').isel(radius=np.arange(0, 30))*2 # convert half hourly to hourly
gpm_polar_mean = get_precip_gpm(cropped_gpm, var_name = 'precip', dates=['27', '28'], hurr_center_ref = True, obs=harvey)['precip'].mean(dim='angle')*2

import matplotlib.pyplot as plt
import numpy as np

# Assuming pre_precip_polar_mean is your xarray DataArray with polar coordinates

# Create the figure and axes
fig, ax = plt.subplots(2, 2, figsize=(10, 9), subplot_kw={'projection': 'polar'},
    gridspec_kw={'bottom': 0.2, 'left':0.09, 'right':0.81, 'hspace':0.3})

cont1 = ax[0, 0].contourf(gpm_polar.angle, gpm_polar.radius, gpm_polar.T, cmap='jet', extend='both', levels=levels)
ax[0, 0].set_theta_zero_location('N')
ax[0, 0].set_title('GPM')

cont2 = ax[1, 0].contourf(pre_precip_polar_mean.angle, pre_precip_polar_mean.radius, pre_precip_polar_mean.T, cmap='jet', extend='both', levels=levels)
ax[1, 0].set_theta_zero_location('N')
ax[1, 0].set_title('LULC 2001')

cont3 = ax[1, 1].contourf(post_precip_polar_mean.angle, post_precip_polar_mean.radius, post_precip_polar_mean.T, cmap='jet', extend='both', levels=levels)
ax[1, 1].set_theta_zero_location('N')
ax[1, 1].set_title('LULC 2017')

ax_reg = fig.add_subplot(2, 2, 2)  
pre_precip_polar_mean.mean(dim='angle').plot(ax=ax_reg, color='b', label='LULC 2001')
post_precip_polar_mean.mean(dim='angle').plot(ax=ax_reg, color='r', label='LULC 2017')
gpm_polar.mean(dim='angle').plot(ax=ax_reg, color='k', label='GPM')
ax_reg.set_xticks(np.arange(0, 300, 60))
ax_reg.legend(fontsize=12)
ax_reg.grid(True)
ax_reg.set_xlim((0, 240))
ax_reg.set_xlabel('Radius from the storm center (km)')
ax_reg.set_ylabel('Precipitation Rate (mm/hr)')


cax = fig.add_axes([0.87, 0.18, 0.03, 0.7])
cbar = fig.colorbar(cont3, cax=cax, shrink=0.75)
cbar.ax.set_ylabel('Precipitation Rate (mm/hr)')

[ax.set_yticks(np.arange(0, 260, 50)) for ax in (ax[0, 0], ax[1, 0], ax[1, 1])]
[ax.tick_params(axis='y', labelsize=12) for ax in (ax[0, 0], ax[1, 0], ax[1, 1])]


ax[0, 0].set_xticklabels([])
ax[1, 0].set_xticklabels([])
ax[1, 1].set_xticklabels([])
ax[0, 1].set_xticklabels([])

plt.tight_layout()
plt.savefig('../figures_draft01/extended_figure02_polar_plot.jpeg', dpi=400, bbox_inches='tight')


fig, axs = plt.subplots(1, 3, figsize=(19, 7), sharex=True, sharey=True,
    gridspec_kw={'bottom': 0.2, 'left': 0.12, 'right': 0.85, 'hspace': 0.3})

im2 = gpm_polar_mean.plot(ax=axs[0], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[0].set_title('GPM Precipitation ')
axs[0].set_xlabel('Distance from hurricane center (km)')

im1 = plot_pre.plot(ax=axs[1], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[1].set_title('Precipitation (LULC 2001)')
axs[1].set_xlabel('Distance from hurricane center (km)')

im2 = plot_post.plot(ax=axs[2], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[2].set_title('Precipitation (LULC 2017)')
axs[2].set_xlabel('Distance from hurricane center (km)')


cax = fig.add_axes([0.88, 0.25, 0.012, 0.6])  # Modify the values as needed
cbar = fig.colorbar(im2, cax=cax)
cbar.ax.set_ylabel('Precipitation Rate(mm/hr)')
plt.tight_layout()
plt.savefig('../figures_draft01/extended_figure02_polar_temporal.jpeg', dpi=400, bbox_inches='tight')



