import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import cartopy.crs as ccrs
from coast import plot_coast
import tropycal.tracks as tracks
import matplotlib as mpl

plt.rcParams.update({"font.size": 14, "font.weight": "bold", "savefig.dpi": 300})

basin = tracks.TrackDataset(basin='north_atlantic',source='hurdat',include_btk=False)

hurdat = basin.get_storm(('harvey',2017))

harvey_2001 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/LULC_2001/ARWpost/Harvey_2001_d02.nc"
harvey_2020 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey/LULC_2020/ARWpost/Harvey_2020_d02.nc"


A = xr.open_dataset(harvey_2001)

loc_ids = [
    np.where(
        A["slp"].isel(time=time_ids).squeeze()
        == A["slp"].isel(time=time_ids).squeeze().min()
    )
    for time_ids in range(30)
]

track_lon_2001 = A["lon"][np.array(loc_ids).squeeze()[:, 1]] - 360
track_lat_2001 = A["lat"][np.array(loc_ids).squeeze()[:, 0]]


A = xr.open_dataset(harvey_2020)

loc_ids = [
    np.where(
        A["slp"].isel(time=time_ids).squeeze()
        == A["slp"].isel(time=time_ids).squeeze().min()
    )
    for time_ids in range(30)
]

track_lon_2020 = A["lon"][np.array(loc_ids).squeeze()[:, 1]] - 360
track_lat_2020 = A["lat"][np.array(loc_ids).squeeze()[:, 0]]



cmap = plt.cm.jet
bounds = [0, 32, 64, 83, 96, 113, 137, 150]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# PHYSICS
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
plot_coast(ax)
ax.plot(hurdat['lon'], hurdat['lat'], "k", label="HURDAT ")
ax.plot(track_lon_2001, track_lat_2001, "b", label="LULC 2001")
ax.plot(track_lon_2020, track_lat_2020, "r", label="LULC 2020")
scatter = ax.scatter(hurdat['lon'], hurdat['lat'], c=hurdat['vmax'], cmap=cmap, norm=norm)
plt.legend()
plt.colorbar(scatter)
ax.set_ylim(15, 40)
ax.set_xlim(-110, -85)
plt.savefig('../figures/Track.jpeg')
plt.show()
