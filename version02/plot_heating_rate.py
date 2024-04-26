import matplotlib.pyplot as plt
import glob
import xarray as xr
import numpy as np

plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "font.family": "monospace"}
)



pre = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*nc'))
post = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*nc'))

plot = True

pre_vals, post_vals = [], []

for index in range(len(pre)):
 prefix = post[index].split('/')[-1][:24]
 print(prefix)
 pre_vals.append(xr.open_dataset(pre[index])['H_DIABATIC_interp'].sel(radius=slice(100, 150)).mean('radius').sel(level=slice(1000, 800))*1000)
 post_vals.append(xr.open_dataset(post[index])['H_DIABATIC_interp'].sel(radius=slice(100, 150)).mean('radius').sel(level=slice(1000, 800))*1000)
 if plot:
  fig, axs = plt.subplots(1, 1, figsize=(6, 7))
  axs.plot(pre_vals[-1], pre_vals[-1].level, 'b', marker='o', label='LULC 2001')
  axs.plot(post_vals[-1], post_vals[-1].level, 'r', marker='o', label='LULC 2017')
  axs.set_xlabel(r'Heating Rate X $10^{-3}$ (K/s) ')
  axs.invert_yaxis()
  axs.set_ylabel('Pressure (hPa)')
  plt.legend()
  plt.grid(True)
  plt.tight_layout()
  plt.title(prefix) 
  plt.savefig(f'../../figures_draft02/heating_rate/{prefix}.jpeg', dpi=400, bbox_inches='tight') 
  plt.close()


fig, axs = plt.subplots(1, 1, figsize=(12, 6))
xr.concat(pre_vals, dim='new').T.plot.contourf(ax=axs, levels=np.arange(-2.4, 2.6, 0.2))
plt.gca().invert_yaxis()
plt.title('LULC 2001')

fig, axs = plt.subplots(1, 1, figsize=(12, 6))
xr.concat(post_vals, dim='new').T.plot.contourf(ax=axs, levels=np.arange(-2.4, 2.6, 0.2))
plt.gca().invert_yaxis()
plt.title('LULC 2017')


fig, axs = plt.subplots(1, 1, figsize=(12, 6))
(xr.concat(post_vals, dim='new') -xr.concat(pre_vals, dim='new')).T.plot.contourf(ax=axs, levels=np.arange(-0.5, 0.6, 0.1))
plt.gca().invert_yaxis()
plt.title('LULC 2017')


plt.show()






pre_vals, post_vals = [], []

for index in range(len(pre)):
 prefix = post[index].split('/')[-1][:24]
 print(prefix)
 pre_vals.append(xr.open_dataset(pre[index])['tangential'].sel(radius=slice(0, 100)).mean('radius').sel(level=slice(1000, 700)))
 post_vals.append(xr.open_dataset(post[index])['tangential'].sel(radius=slice(0, 100)).mean('radius').sel(level=slice(1000, 700)))


pre_timeseries = xr.concat(pre_vals, dim='new').mean(dim='level')
post_timeseries = xr.concat(post_vals, dim='new').mean(dim='level')

fig, axs = plt.subplots(1, 1, figsize=(8, 5))
pre_timeseries.plot(ax=axs, color='b', label='LULC 2001')
post_timeseries.plot(ax=axs, color='r', label='LULC 2017')
plt.legend()
plt.show()




