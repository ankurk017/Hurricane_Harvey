import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


R1_file = '/nas/rstor/akumar/CM1/cm1r20.2/run/cm10ut/cm1out_s_roughness01.nc'
R2_file = '/nas/rstor/akumar/CM1/cm1r20.2/run/cm10ut/cm1out_s_roughness02.nc'

R1 = xr.open_dataset(R1_file)
R2 = xr.open_dataset(R2_file)

fig, axs = plt.subplots(2, 1, figsize=(8, 12))

cont1 = axs[0].contourf(R1['lon'], np.arange(0, 151, 1), R1['rain'].squeeze().values, levels=np.arange(100, 1100, 100), cmap='jet', extend='both')
cont2 = axs[1].contourf(R2['lon'], np.arange(0, R2['rain'].shape[0], 1), R2['rain'].squeeze().values, levels=np.arange(100, 1100, 100), cmap='jet', extend='both')
[ax.set_xlim(0, 60) for ax in axs]
[ax.set_ylim(0, 150) for ax in axs]

[ax.set_xlabel('Radius (km)') for ax in axs]
[ax.set_ylabel('Forecast Hours (hr)') for ax in axs]
[ax.set_yticks(np.arange(0, 150, 24)) for ax in axs]

cbar = plt.colorbar(cont2, ax=axs.ravel())
cbar.ax.set_ylabel('Accumulated Rainfall (mm)', rotation=90, labelpad=15)

fig, axs = plt.subplots(2, 1, figsize=(8, 12))

cont1 = axs[0].contourf(R1['lon'], np.arange(0, 151, 1), np.sqrt(R1['u10'].squeeze().values**2+R1['v10'].squeeze().values**2), levels=np.arange(0, 40, 4), cmap='jet', extend='both')
cont2 = axs[1].contourf(R2['lon'], np.arange(0, R2['rain'].shape[0], 1), np.sqrt(R2['u10'].squeeze().values**2+R2['v10'].squeeze().values**2), levels=np.arange(0, 40, 4), cmap='jet', extend='both')
[ax.set_xlim(0, 60) for ax in axs]
[ax.set_ylim(0, 150) for ax in axs]

[ax.set_xlabel('Radius (km)') for ax in axs]
[ax.set_ylabel('Forecast Hours (hr)') for ax in axs]
[ax.set_yticks(np.arange(0, 150, 24)) for ax in axs]

cbar = plt.colorbar(cont2, ax=axs.ravel())
cbar.ax.set_ylabel('10 m Wind Speed (m/s)', rotation=90, labelpad=15)

fig, axs = plt.subplots(2, 1, figsize=(8, 12))

cont1 = axs[0].contourf(R1['lon'], np.arange(0, 151, 1), R1['ust'].squeeze().values, levels=np.arange(0.5, 5, 0.5), cmap='gist_ncar', extend='both')
cont2 = axs[1].contourf(R2['lon'], np.arange(0, R2['rain'].shape[0], 1), R2['ust'].squeeze().values, levels=np.arange(0.5, 5, 0.5), cmap='gist_ncar', extend='both')
[ax.set_xlim(0, 60) for ax in axs]
[ax.set_ylim(0, 150) for ax in axs]

[ax.set_xlabel('Radius (km)') for ax in axs]
[ax.set_ylabel('Forecast Hours (hr)') for ax in axs]
[ax.set_yticks(np.arange(0, 150, 24)) for ax in axs]

cbar = plt.colorbar(cont2, ax=axs.ravel())
cbar.ax.set_ylabel('Friction Velocity (m/s)', rotation=90, labelpad=15)

plt.show()



