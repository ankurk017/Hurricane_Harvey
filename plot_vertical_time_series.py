import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset

from wrf import to_np, getvar, CoordPair, vertcross
from src.wrf_src import wrf_assign_coords
import glob
import progressbar
import xarray as xr

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


case = "pre"
home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/"
pre_files = sorted(
    glob.glob(
        home_2512
        + f"/WRF_FNL_2612/pre_UCM/WRF_2dom/test/em_real/wrfout_d01_2017-*"
    )
)[::3]

case = "post"
post_files = sorted(
    glob.glob(
        home_2512
        + f"/WRF_FNL_2612/post_UCM/WRF_2dom/test/em_real/wrfout_d01_2017-*"
    )
)[::3]


location = (-95.499,  29.74)  # Houston location
box = 0.2
pre_w_loc = []
post_w_loc = []

for fileid in progressbar.progressbar(range(len(post_files))):

    pre_ncfile = Dataset(pre_files[fileid])
    pre_wspd = wrf_assign_coords(getvar(pre_ncfile, "wa", units="m/s"))
    post_ncfile = Dataset(post_files[fileid])
    post_wspd = wrf_assign_coords(getvar(post_ncfile, "wa", units="m/s"))

    lat_id = np.where(
        np.logical_and(
            pre_wspd["south_north"].values > location[1] - box,
            pre_wspd["south_north"].values <= location[1] + box,
        )
    )[0]
    lon_id = np.where(
        np.logical_and(
            pre_wspd["west_east"].values > location[0] - box,
            pre_wspd["west_east"].values <= location[0] + box,
        )
    )[0]


    pre_w_loc.append(pre_wspd.isel(south_north=lat_id, west_east=lon_id).mean(dim='south_north').mean(dim='west_east'))
    post_w_loc.append(post_wspd.isel(south_north=lat_id, west_east=lon_id).mean(dim='south_north').mean(dim='west_east'))

pre_vert_vel = xr.concat(pre_w_loc, dim="Time").T
post_vert_vel = xr.concat(post_w_loc, dim="Time").T
vert_cel_diff = post_vert_vel-pre_vert_vel

z_loc = wrf_assign_coords(getvar(pre_ncfile, "z")).isel(south_north=lat_id, west_east=lon_id).mean(dim='south_north').mean(dim='west_east')/1000

plt.figure(figsize=(10, 5))
plot = vert_cel_diff.assign_coords({"bottom_top": z_loc.values}).plot()
plot.colorbar.set_label('Vertical Velocity (m/s)')
plt.ylabel('Height (km)')
plt.title('Post - Pre Vertical Velocity')
plt.tight_layout()
plt.show()



