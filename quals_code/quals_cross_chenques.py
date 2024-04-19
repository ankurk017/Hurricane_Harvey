import progressbar
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
import glob


wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/"
wrf_post = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[::3][:18]

wrf_runs = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
wrf_pre = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[::3][:18]

start = (22., -92.5)
end = (28, -96.9)
pres_levels = np.arange(1000, 100, -10)

uv_cs = []
for files in progressbar.progressbar(wrf_pre):
    ncfile = Dataset(files)
    slp = getvar(ncfile, "wspd_wdir10").sel(wspd_wdir='wspd')
    slp = getvar(ncfile, "slp")
    p = getvar(ncfile, "pressure")
    z = getvar(ncfile, "z", units="dm")
    uv = slp
    lats, lons = latlon_coords(slp)
    cart_proj = get_cartopy(slp)
    uv_xr = wrf_assign_coords(uv.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    ))
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})
    uv_cs.append(cross_section(uv_xr, start, end))

uv_cs_post = []
for files in progressbar.progressbar(wrf_post):
    ncfile = Dataset(files)
    slp = getvar(ncfile, "wspd_wdir10").sel(wspd_wdir='wspd')
    slp = getvar(ncfile, "slp")
    p = getvar(ncfile, "pressure")
    z = getvar(ncfile, "z", units="dm")
    uv = slp
    lats, lons = latlon_coords(slp)
    cart_proj = get_cartopy(slp)
    uv_xr = wrf_assign_coords(uv.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    ))
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})
    uv_cs_post.append(cross_section(uv_xr, start, end))

uvcs = xr.concat(uv_cs_post, dim='Time') - xr.concat(uv_cs, dim='Time')

dis = np.sqrt((uvcs['x'].values-(-96.9))**2 +
              (uvcs['y'].values-(28))**2)*111.11
uvcs_dist = uvcs.assign_coords({'index': dis})
