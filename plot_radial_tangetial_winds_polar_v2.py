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

wrf_file = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_Quals/WRF_cropland/WRF/wrfout_d02_2017-08-26_00:00:00"


pres_levels = np.arange(1000, 100, -10)

ncfile = Dataset(wrf_file)

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

uv_xr = uv_xr.assign_coords(lat=uv_xr['lat'].values-cen_loc[1])
uv_xr = uv_xr.assign_coords(lon=uv_xr['lon'].values-cen_loc[0])

r = np.arange(0, 2, 0.1)
ang = np.arange(0 ,361, 1)*np.pi/180
r_m, ang_m = np.meshgrid(r, ang)

x = r_m*np.cos(ang_m)
y = r_m*np.sin(ang_m)


start = (27.61, -96.75)
end = (27.61, -98.31)


uv_cs = cross_section(uv_xr, start, end).set_coords(('y', 'x'))
uv_cs['radial'], uv_cs['tangential'] = metpy.calc.cross_section_components(
    uv_cs['ua_interp'], uv_cs['va_interp'])


fig, axs = plt.subplots(1, 2, figsize=(12, 5))
uv_cs['ua_interp'].plot(ax=axs[0], yincrease=False)
uv_cs['va_interp'].plot(ax=axs[1], yincrease=False)
[ax.set_yscale("log") for ax in axs]
plt.tight_layout()

fig, axs = plt.subplots(1, 2, figsize=(12, 5))
uv_cs['radial'].plot(ax=axs[0], yincrease=False)
uv_cs['tangential'].plot(ax=axs[1], yincrease=False)
[ax.set_yscale("log") for ax in axs]
plt.tight_layout()

plt.show()


