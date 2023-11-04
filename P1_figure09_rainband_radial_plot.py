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



def convert_polar1(var_cropped, margin=5, resolution=0.09):
    r = np.arange(0, margin, resolution)
    ang = np.arange(0, 361, 1) * np.pi / 180
    r_mesh, ang_mesh = np.meshgrid(r, ang)

    x_polar = r_mesh * np.cos(ang_mesh)
    y_polar = r_mesh * np.sin(ang_mesh)

    xlong_values = var_cropped["south_north"]
    xlat_values = var_cropped["west_east"]

    pcp_values = var_cropped.values

    pcp_interp = RegularGridInterpolator((xlong_values, xlat_values), pcp_values, method='linear', bounds_error=False, fill_value=np.nan)

    x_polar_flat = x_polar.flatten()
    y_polar_flat = y_polar.flatten()

    pcp_polar_flat = pcp_interp((x_polar_flat, y_polar_flat))
    pcp_polar = pcp_polar_flat.reshape(x_polar.shape)


    return r, ang, pcp_polar

def get_precip(wrf_runs, var_name = 'precip', dates=['27', '28'], hurr_center_ref = True, ref_loc = (-97.2, 29.2)):

    wrfoutfile = [sorted(glob.glob(wrf_runs + f"wrfout_d02*-{date}_*")) for date in dates]
    wrfoutfile = list(itertools.chain.from_iterable(wrfoutfile))

    basin = tracks.TrackDataset(
        basin="north_atlantic", source="hurdat", include_btk=False
    )

    harvey = basin.get_storm(("harvey", 2017))
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


pre_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
post_files = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
variable = 'precip'

pre_precip = get_precip(pre_files, dates=['27', '28'], var_name = variable)[variable]
post_precip = get_precip(post_files, dates=['27', '28'], var_name = variable)[variable]


pre_precip_polar_mean = pre_precip.mean(dim='time').isel(radius=np.arange(0, 25))
post_precip_polar_mean =  post_precip.mean(dim='time').isel(radius=np.arange(0, 25))

minmax = find_common_min_max((pre_precip_polar_mean.values, post_precip_polar_mean.values))
levels = np.linspace(minmax[0], minmax[1], 15)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': 'polar'},
    gridspec_kw={'bottom': 0.2, 'left':0.09, 'right':0.81, 'hspace':0.3})

#cont1 = ax1.contourf(ang, r[:radius_bound] * 111.11, pre_precip_polar_mean[:, :radius_bound].T, cmap='jet', extend='both', levels=levels)
cont1 = ax1.contourf(pre_precip_polar_mean.angle, pre_precip_polar_mean.radius, pre_precip_polar_mean.T, cmap='jet', extend='both', levels=levels)
ax1.set_theta_zero_location('N')
ax1.set_title('LULC 2001')

#cont2 = ax2.contourf(ang, r[:radius_bound] * 111.11, post_precip_polar_mean[:, :radius_bound].T, cmap='jet', extend='both', levels=levels)
cont2 = ax2.contourf(post_precip_polar_mean.angle, post_precip_polar_mean.radius, post_precip_polar_mean.T, cmap='jet', extend='both', levels=levels)
ax2.set_theta_zero_location('N')
ax2.set_title('LULC 2017')

cax = fig.add_axes([0.87, 0.1, 0.03, 0.8])
cbar = fig.colorbar(cont2, cax=cax, shrink=0.75)
cbar.ax.set_ylabel('Precipitation (mm/hr)')

# Remove radial grid labels
#ax1.set_yticklabels([])
#ax2.set_yticklabels([])

ax1.set_xticklabels([])
ax2.set_xticklabels([])

plt.tight_layout()
plt.savefig(f'../figures_paper/rainband/pcp_polar.jpeg', dpi=300)


pre_precip = get_precip(pre_files, dates=['27', '28'], var_name = variable, hurr_center_ref=True, ref_loc=(-96.8, 29.2))[variable]
post_precip = get_precip(post_files, dates=['27', '28'], var_name = variable, hurr_center_ref=True, ref_loc=(-96.8, 29.2))[variable]

plot_pre = pre_precip.mean(dim='angle')
plot_post = post_precip.mean(dim='angle')

minmax = find_common_min_max((plot_pre.values, plot_post.values))
levels = np.linspace(minmax[0], minmax[1]*0.8, 15)

fig, axs = plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True,
    gridspec_kw={'bottom': 0.2, 'left': 0.12, 'right': 0.85, 'hspace': 0.3})

im1 = plot_pre.plot(ax=axs[0], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[0].set_title('Precipitation (LULC 2001)')
axs[0].set_xlabel('Distance from hurricane center (km)')

im2 = plot_post.plot(ax=axs[1], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[1].set_title('Precipitation (LULC 2017)')
axs[1].set_xlabel('Distance from hurricane center (km)')
cax = fig.add_axes([0.88, 0.25, 0.02, 0.6])  # Modify the values as needed
cbar = fig.colorbar(cont2, cax=cax)
cbar.ax.set_ylabel('Precipitation (mm/hr)')
plt.savefig(f'../figures_paper/rainband/Rainfall_TS_radial_mean.jpeg', dpi=300)

plot_pre = pre_precip.sel(angle=slice(0, 2)).mean(dim='angle')
plot_post = post_precip.sel(angle=slice(0, 2)).mean(dim='angle')

minmax = find_common_min_max((plot_pre.values, plot_post.values))
levels = np.linspace(minmax[0], minmax[1]*0.8, 15)

fig, axs = plt.subplots(1, 2, figsize=(14, 7), sharex=True, sharey=True,
    gridspec_kw={'bottom': 0.2, 'left': 0.12, 'right': 0.85, 'hspace': 0.3})

im1 = plot_pre.plot(ax=axs[0], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[0].set_title('Precipitation (LULC 2001)')
axs[0].set_xlabel('Distance from hurricane center (km)')

im2 = plot_post.plot(ax=axs[1], cmap='jet', levels=levels, extend='both', add_colorbar=False)
axs[1].set_title('Precipitation (LULC 2017)')
axs[1].set_xlabel('Distance from hurricane center (km)')
cax = fig.add_axes([0.88, 0.25, 0.02, 0.6])  # Modify the values as needed
cbar = fig.colorbar(cont2, cax=cax)
cbar.ax.set_ylabel('Precipitation (mm/hr)')
plt.savefig(f'../figures_paper/rainband/Rainfall_TS_radial_sliced.jpeg', dpi=300)



time_slice = ['2017-08-27', '2017-08-28']


fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), sharex=True, sharey=True)

for index in range(2):
 pre_mean = pre_precip.sel(time=time_slice[index]).sel(angle=slice(0, 2)).mean(dim=['time', 'angle'])
 post_mean = post_precip.sel(time=time_slice[index]).sel(angle=slice(0, 2)).mean(dim=['time', 'angle'])
 pre_mean.plot(label='LULC 2001', color='b',ax=ax.ravel()[index])
 post_mean.plot(label='LULC 2017', color='r', ax=ax.ravel()[index])

 ax.ravel()[index].set_title(f'Radial mean of precipitaion for {time_slice[index]}')
 ax.ravel()[index].set_xlabel('Radius from storm center (km)')
 ax.ravel()[index].set_ylabel('Mean Precipitation (mm/hr)')
 ax.ravel()[index].grid()
ax[1].legend()

plt.tight_layout()
plt.savefig(f'../figures_paper/rainband/radial_pcp.jpeg', dpi=300)
#plt.show()




"""
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='polar')
cont = ax.contourf(ang, r[:15]*111.11, post_precip_polar_mean[:,
                   :15].T - pre_precip_polar_mean[:, :15].T, cmap='bwr')
ax.set_theta_zero_location('N')
cbar = plt.colorbar(cont, ax=ax)
cbar.ax.set_ylabel('Precipitation Error (mm/hr)')
plt.tight_layout()


plt.figure(figsize=(7, 4))
plt.plot(r * 111.111, np.nanmean(pre_precip_polar_mean.T, axis=1),
         'r', label='LULC 2001')
plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T, axis=1),
         'b', label='LULC 2017')
plt.xlabel("Radius from the storm center (km)")
plt.ylabel("Precipation")
plt.tight_layout()
plt.legend()
plt.xlim((0, np.round(r.max() * 111.11)))


plt.figure(figsize=(7, 4))
plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T -
         pre_precip_polar_mean.T, axis=1), 'r', label='LULC 2001')
# plt.plot(r * 111.111, np.nanmean(post_precip_polar_mean.T, axis=1), 'b', label='LULC 2017')
plt.xlabel("Radius from the storm center (km)")
plt.ylabel("Precipation")
plt.tight_layout()
plt.legend()
plt.xlim((0, np.round(r.max() * 111.11)))
plt.show()
"""
