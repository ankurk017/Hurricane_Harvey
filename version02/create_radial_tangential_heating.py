import matplotlib.pyplot as plt
import glob
from netCDF4 import Dataset
from wrf import getvar, interplevel
import numpy as np
from src.wrf_src import wrf_assign_coords
import metpy
from metpy.interpolate import cross_section
import xarray as xr
import progressbar


from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

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

# from self_utils import ahps, coast
import src.coast as coast
import src.ahps as ahps
import numpy as np
from scipy.interpolate import griddata
from matplotlib import gridspec
from src.wrf_src import wrf_assign_coords
import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
import xarray as xr
import metpy
from metpy.interpolate import cross_section
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import progressbar


plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "font.family": "monospace"}
)


def get_radial_tangential_wind(infile, pres_levels=np.concatenate((np.arange(1000, 790, -10), np.arange(750, 100, -50))), radial_extent=5
                               ):

    ncfile = Dataset(infile)

    slp = wrf_assign_coords(getvar(ncfile, "slp"))
    cen_lat = np.unique(slp.south_north[np.where(slp == slp.min())[0]].values)
    cen_lon = np.unique(slp.west_east[np.where(slp == slp.min())[1]].values)
    cen_loc = (cen_lon[0], cen_lat[0])

    p = getvar(ncfile, "pressure")

    ua = getvar(ncfile, "ua", units="m/s")
    va = getvar(ncfile, "va", units="m/s")
    wa = getvar(ncfile, "wa", units="m/s")
    hr = getvar(ncfile, "H_DIABATIC")

    ua_pres = xr.concat([interplevel(ua, p, pres)
                        for pres in pres_levels], dim="level")
    va_pres = xr.concat([interplevel(va, p, pres)
                        for pres in pres_levels], dim="level")
    wa_pres = xr.concat([interplevel(wa, p, pres)
                        for pres in pres_levels], dim="level")
    hr_pres = xr.concat([interplevel(hr, p, pres)
                        for pres in pres_levels], dim="level")

    uv = xr.merge((ua_pres, va_pres, wa_pres, hr_pres))

    uv_xr = wrf_assign_coords(
        uv.metpy.assign_crs(
            grid_mapping_name="latitude_longitude", earth_radius=6371229.0
        )
    )
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})

    uv_xr["y"] = uv_xr["y"].values - cen_loc[1]
    uv_xr["x"] = uv_xr["x"].values - cen_loc[0]
    uv_xr = uv_xr.drop('Time')
    start = (0, 0)
    theta = np.linspace(0, np.pi / 4, 8)
    # theta = np.linspace(0, 2*np.pi, 16)

    x = radial_extent * np.cos(theta)
    y = radial_extent * np.sin(theta)

    uv_rad_tang = []
    for index in range(len(x)):
        uv_cs = cross_section(
            uv_xr, start, (y[index], x[index])).set_coords(("y", "x"))
        uv_cs["radial"], uv_cs["tangential"] = metpy.calc.cross_section_components(
            uv_cs["ua_interp"], uv_cs["va_interp"]
        )
        uv_cs["index"] = (
            (np.sqrt(uv_cs["y"] ** 2 + uv_cs["x"] ** 2)) * 111.11).values
        uv_cs = uv_cs.rename({'index': 'radius'})
        uv_cs = uv_cs.interp(radius=np.arange(0, 500, 5))
        uv_rad_tang.append(uv_cs)

    uv_radial_tangential_wind = xr.concat(
        uv_rad_tang, dim="new").mean(dim="new")
    return uv_radial_tangential_wind  # .drop('metpy_crs')


wrf_files_pre = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*00'))[:85]
wrf_files_post = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-*00'))[:85]


for index in progressbar.progressbar(range(len(wrf_files_post))):
    prefix = wrf_files_post[index].split('/')[-1][11:24]
    # print(prefix)
    levels = np.concatenate(
        (np.arange(1000, 790, -10), np.arange(750, 100, -50)))

    get_radial_tangential_wind(wrf_files_pre[index], ).drop('metpy_crs') .to_netcdf(
        wrf_files_pre[index][:-6]+'_radial_tangential_heating.nc')
    get_radial_tangential_wind(wrf_files_post[index], ).drop('metpy_crs').to_netcdf(
        wrf_files_post[index][:-6]+'_radial_tangential_heating.nc')
