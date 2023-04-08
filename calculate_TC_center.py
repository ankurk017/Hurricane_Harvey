import numpy as np
import xarray as xr
import glob
import matplotlib.pyplot as plt
wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/'
wrfoutfile_pre = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

from netCDF4 import Dataset
from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim)


import tropycal.tracks as tracks

from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

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
from src.wrf_src import wrf_assign_coords
import metpy
from metpy.interpolate import cross_section
import numpy as np
import xarray as xr

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})



f = wrfoutfile_pre[12]

ncfile = Dataset(f)

basin = tracks.TrackDataset(basin='north_atlantic',
                            source='hurdat', include_btk=False)
harvey = basin.get_storm(('harvey', 2017))
lat = lats.values[:,0]

# Extract the pressure, geopotential height, and wind variables
p = getvar(ncfile, "pressure")
z = getvar(ncfile, "z", units="dm")
ua = getvar(ncfile, "ua", units="kt")
va = getvar(ncfile, "va", units="kt")
wspd = getvar(ncfile, "wspd_wdir", units="kts")[0,:]

uv10 = getvar(ncfile, 'uvmet10')
u10 = uv10.sel(u_v='u').values
v10 = uv10.sel(u_v='v').values
slp = getvar(ncfile, 'slp')
lats, lons = latlon_coords(slp)
slp = slp.values
lat = lats.values[:,0]
lon = lons.values[0,:]

obs_loc = [-93.1, 23.7 ]
min_pressure_index = np.unravel_index(np.argmin(slp), slp.shape)
min_pressure_location = (lon[min_pressure_index[1]], lat[min_pressure_index[0]])
print("Location of minimum pressure:", min_pressure_location)
pressure = slp[79:179, 79:179]
x = lon[79:179]
y = lat[79:179]
x_grid, y_grid = np.meshgrid(x, y) 
x_centroid = np.sum(pressure * x_grid) / np.sum(pressure)
y_centroid = np.sum(pressure * y_grid) / np.sum(pressure)
pressure_centroid = (x_centroid, y_centroid)
print("Pressure centroid:", pressure_centroid)

plt.figure(figsize=(7, 6))
#plt.contourf(lon, lat, slp)
plt.imshow(pressure, extent=[x.min(), x.max(), y.min(), y.max()])
plt.plot(obs_loc[0], obs_loc[1], 'k+')
plt.scatter(min_pressure_location[0], min_pressure_location[1], color='r')
plt.scatter(pressure_centroid[0], pressure_centroid[1], color='b')




start = (23.6, -93.9)
end = (23.6, -92.4)
pres_levels = np.arange(1000, 100, -10)

ncfile = Dataset(f)

slp = getvar(ncfile, "slp")
p = getvar(ncfile, "pressure")
z = getvar(ncfile, "z", units="dm")

ua = getvar(ncfile, "ua", units="kt")
va = getvar(ncfile, "va", units="kt")

ua_pres = xr.concat([interplevel(ua, p, pres)
                    for pres in pres_levels], dim='level')
va_pres = xr.concat([interplevel(va, p, pres)
                    for pres in pres_levels], dim='level')

uv = xr.merge((ua_pres, va_pres))


lats, lons = latlon_coords(slp)

cart_proj = get_cartopy(slp)


uv_xr = wrf_assign_coords(uv.metpy.assign_crs(
    grid_mapping_name="latitude_longitude", earth_radius=6371229.0
))
uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})

uv_cs = cross_section(uv_xr, start, end).set_coords(('y', 'x'))
uv_cs['radial'], uv_cs['tangential'] = metpy.calc.cross_section_components(
    uv_cs['ua_interp'], uv_cs['va_interp'])


fig, axs = plt.subplots(1, 2, figsize=(12, 5))
uv_cs['ua_interp'].plot(ax=axs[0], yincrease=False)
uv_cs['va_interp'].plot(ax=axs[1], yincrease=False)
[ax.set_yscale("log") for ax in axs]
plt.tight_layout()


#fig, axs = plt.subplots(1, 2, figsize=(12, 5))
#uv_cs['radial'].plot(ax=axs[0], yincrease=False)
#uv_cs['tangential'].plot(ax=axs[1], yincrease=False)
#[ax.set_yscale("log") for ax in axs]
#plt.tight_layout()

ws = np.sqrt(uv_cs['ua_interp']**2+uv_cs['va_interp']**2).values
fig, axs = plt.subplots(1, 1, figsize=(8, 7))
plt.contourf(uv_cs['x'].values, uv_cs['level'].values, ws, cmap='jet')
plt.scatter(pressure_centroid[0], 1000, color='b')
plt.scatter(min_pressure_location[0], 1000, color='r')
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()

ws = np.sqrt(uv_cs['ua_interp']**2).values
fig, axs = plt.subplots(1, 1, figsize=(8, 7))
plt.contourf(uv_cs['x'].values, uv_cs['level'].values, ws, cmap='bwr')
plt.scatter(pressure_centroid[0], 1000, color='b')
plt.scatter(min_pressure_location[0], 1000, color='r')
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()







