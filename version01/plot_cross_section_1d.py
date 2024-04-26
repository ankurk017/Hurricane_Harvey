from wrf import (
    getvar,
)
from src.wrf_src import wrf_assign_coords
from metpy.interpolate import cross_section
import numpy as np
import xarray as xr
import glob
from wrf import CoordPair

import progressbar
from netCDF4 import Dataset

# cross section 2
start_point = CoordPair(lat=29.3, lon=-96.5)
end_point = CoordPair(lat=30.5, lon=-94.5)


def to_metpy(uv):
    uv_xr = wrf_assign_coords(
        uv.metpy.assign_crs(
            grid_mapping_name="latitude_longitude", earth_radius=6371229.0
        )
    )
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})
    return uv_xr


home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"

wrf_pre = sorted(glob.glob(home_2512 + f"/pre_UCM/WRF/test/em_real/wrfout_d02_2017-*"))[
    :36
]
wrf_post = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF/test/em_real/wrfout_d02_2017-*")
)[:36]


start = (start_point.lat, start_point.lon)
end = (end_point.lat, end_point.lon)
pres_levels = np.arange(1000, 100, -10)

uv = []
for files in progressbar.progressbar(wrf_pre):
    ncfile = Dataset(files)
    uv.append(to_metpy(getvar(ncfile, "wspd_wdir10").sel(wspd_wdir="wspd")))
#    uv.append(to_metpy(getvar(ncfile, "slp")))
uv_cs_pre = cross_section(xr.concat(uv, dim="Time"), start, end)

uv = []
for files in progressbar.progressbar(wrf_post):
    ncfile = Dataset(files)
    uv.append(to_metpy(getvar(ncfile, "wspd_wdir10").sel(wspd_wdir="wspd")))
#    uv.append(to_metpy(getvar(ncfile, "slp")))
uv_cs_post = cross_section(xr.concat(uv, dim="Time"), start, end)


xticks = [str(np.round(vals[0].values, 2))+ ', '+str(np.round(vals[1].values, 2)) for vals in zip(uv_cs.coords['x'], uv_cs.coords['y'])]

plt.figure(figsize=(10, 7))
(uv_cs_post - uv_cs_pre).plot.contourf(cmap='bwr')
#plt.contourf(np.arange(0, 258), uv_cs.Time, uv_cs)
plt.show()





