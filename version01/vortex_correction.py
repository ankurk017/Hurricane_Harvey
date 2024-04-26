import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from src.wrf_src import wrf_assign_coords
from src.coast import plot_coast
import cartopy.crs as ccrs

input_file = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/wrfinput_d02_raw'

A = xr.open_dataset(input_file)


#replace_west_east_id = abs(slp.west_east-(95.8)).argmin()
#replace_south_north_id = abs(slp.south_north-26.3).argmin()

#A['PSFC'].values[0, replace_south_north_id, replace_west_east_id] = 94900

slp = A['PSFC']*0.97
u10 = A['U10']*1.5
v10 = A['V10']*1.5

A['U10'] = u10
A['V10'] = v10
A['PSFC'] = slp

A.to_netcdf('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/wrfinput_d02_vortex')

raw = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/wrfinput_d01_raw'
vortex = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/wrfinput_d01_vortex'
wrfda = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/wrfinput_d01_wrfda'

f_raw = xr.open_dataset(raw)
f_vortex = xr.open_dataset(vortex)
f_wrfda = xr.open_dataset(wrfda)

slice_lonlat = {'south_north':slice(25, 29), 'west_east':slice(-98, -94)}
var = 'PSFC'
fig, ax = plt.subplots(1, 3, figsize=(18, 5), subplot_kw={'projection': ccrs.PlateCarree()})

wrf_assign_coords(f_raw[var].squeeze()).sel(slice_lonlat).plot(ax=ax[0], cmap='jet', cbar_kwargs={'shrink':0.8})
wrf_assign_coords(f_vortex[var].squeeze()).sel(slice_lonlat).plot(ax=ax[1], cmap='jet', cbar_kwargs={'shrink':0.8})
wrf_assign_coords(f_wrfda[var].squeeze()).sel(slice_lonlat).plot(ax=ax[2], cmap='jet', cbar_kwargs={'shrink':0.8})
[plot_coast(ax, houston=True, linewidth=1.5) for ax in ax]


plt.show()









