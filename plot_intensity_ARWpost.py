import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import cartopy.crs as ccrs
import tropycal.tracks as tracks
import matplotlib as mpl

plt.rcParams.update({"font.size": 14, "font.weight": "bold", "savefig.dpi": 300})

basin = tracks.TrackDataset(basin='north_atlantic',source='hurdat',include_btk=False)

hurdat = basin.get_storm(('harvey',2017))

harvey_2001 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_GFS/LULC_2001/ARWpost/Harvey_2001_d02.nc"
harvey_2020 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/WRF_GFS/LULC_2020/ARWpost/Harvey_2020_d02.nc"


A = xr.open_dataset(harvey_2001)

mslp_2001 = np.array([
        A["slp"].isel(time=time_ids).squeeze().min().values
    for time_ids in range(30)
])
ws_2001 = np.array([
        A["ws10"].isel(time=time_ids).squeeze().max().values
    for time_ids in range(30)
])*2

A = xr.open_dataset(harvey_2020)

mslp_2020 = np.array([
        A["slp"].isel(time=time_ids).squeeze().min().values
    for time_ids in range(30)
])
ws_2020 = np.array([
        A["ws10"].isel(time=time_ids).squeeze().max().values
    for time_ids in range(30)
])*2



plt.figure(figsize=(8, 5))
plt.plot(hurdat['date'], hurdat['mslp'], 'k', label='HURDAT')
plt.plot(A['time'][:30], mslp_2001, 'b', label='LULC_2001')
plt.plot(A['time'][:30], mslp_2020, 'r', label='LULC_2020')
plt.legend(fontsize=12)
plt.xticks(rotation = 45) 
plt.tight_layout()
#plt.savefig('../figures/Intensity_MSLP.jpeg')


plt.figure(figsize=(8, 5))
plt.plot(hurdat['date'], hurdat['vmax'], 'k', label='HURDAT')
plt.plot(A['time'][:30], ws_2001, 'b', label='LULC_2001')
plt.plot(A['time'][:30], ws_2020, 'r', label='LULC_2020')
plt.legend(fontsize=12)
plt.xticks(rotation = 45) 
plt.tight_layout()
#plt.savefig('../figures/Intensity_WS.jpeg')
plt.show()


