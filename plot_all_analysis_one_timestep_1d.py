from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np
from scipy.interpolate import griddata
import progressbar
import matplotlib
from src.wrf_src import wrf_assign_coords, crop_region, plot_bb
# matplotlib.use('Agg')
import xarray as xr


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})

lons=[-95.35, -95.82, -95.70, -95.46, -95.37] 
lats=[ 29.75,  29.73,  29.91,  30.01, 29.57]

lons=[-95.35, -95.82, ] 
lats=[ 29.75,  29.73,]

location1 = {"south_north": lats[0], "west_east": lons[0], "box": 0.2}
location2 = {"south_north": lats[1], "west_east": lons[1], "box": 0.2}


home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/"
wrfoutfile_pre = sorted(
    glob.glob(home_2512 + f"/pre_UCM/WRF/test/em_real/wrfout_d02_2017-*")
)[:50]
wrfoutfile_post = sorted(
    glob.glob(home_2512 + f"/post_UCM/WRF/test/em_real/wrfout_d02_2017-*")
)[:50]

index = np.min((len(wrfoutfile_pre), len(wrfoutfile_post)))

prefiles = wrfoutfile_pre[:index]
postfiles = wrfoutfile_post[:index]


years = (
    2001,
    2017,
)

plot_bb(prefiles[0], (location1, location2))



wrf_pcp_pre_windspeed_loc1, wrf_pcp_post_windspeed_loc1, wrf_pcp_pre_windspeed_loc2, wrf_pcp_post_windspeed_loc2 = [], [], [], []

for fileid in progressbar.progressbar(range(len(prefiles) - 1)):
    wrf_pcp_pre_windspeed_loc1.append(crop_region(
        #wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "wspd_wdir10").sel(wspd_wdir='wspd')),
        wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "T2")),
        location=location1,
    ).mean())
    wrf_pcp_post_windspeed_loc1.append(crop_region(
        #wrf_assign_coords(getvar(Dataset(postfiles[fileid + 1]), "wspd_wdir10").sel(wspd_wdir='wspd')),
        wrf_assign_coords(getvar(Dataset(postfiles[fileid + 1]), "T2")),
        location=location1,
    ).mean())

    wrf_pcp_pre_windspeed_loc2.append(crop_region(
        #wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "wspd_wdir10").sel(wspd_wdir='wspd')),
        wrf_assign_coords(getvar(Dataset(prefiles[fileid + 1]), "T2")),
        location=location2,
    ).mean())
    wrf_pcp_post_windspeed_loc2.append(crop_region(
        #wrf_assign_coords(getvar(Dataset(postfiles[fileid + 1]), "wspd_wdir10").sel(wspd_wdir='wspd')),
        wrf_assign_coords(getvar(Dataset(postfiles[fileid + 1]), "T2")),
        location=location2,
    ).mean())

plt.figure(figsize=(8, 5))
xr.concat(wrf_pcp_pre_windspeed_loc1, dim='Time').plot(label='LULC 2001', color='b')
xr.concat(wrf_pcp_post_windspeed_loc1, dim='Time').plot(label='LULC 2017', color='r')
#plt.tight_layout()


#plt.figure(figsize=(8, 5))
xr.concat(wrf_pcp_pre_windspeed_loc2, dim='Time').plot(label='LULC 2001', color='c')
xr.concat(wrf_pcp_post_windspeed_loc2, dim='Time').plot(label='LULC 2017', color='g')
plt.tight_layout()


plt.show()



import src.calculate_lwp as clwp
from metpy.interpolate import cross_section


#start_point = CoordPair(lat=29.55, lon=-96.07)
#end_point = CoordPair(lat=30.27, lon=-94.84)

start = (29.55, -96.07)
end = (30.27, -94.84)


uv_cs_pre = []
for files in progressbar.progressbar(prefiles[3:]):
    ncfile = Dataset(files)
    uv = getvar(ncfile, "wspd_wdir10").sel(wspd_wdir='wspd')
    #clwp_values = clwp.calculate_lwp(filein=files)
    #uv.values = clwp_values
    uv_xr = wrf_assign_coords(uv.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    ))
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})
    uv_xr.Time.values = uv.Time.values
    uv_cs_pre.append(cross_section(uv_xr, start, end))



uv_cs_pst = []
for files in progressbar.progressbar(postfiles[3:]):
    ncfile = Dataset(files)
    uv = getvar(ncfile, "wspd_wdir10").sel(wspd_wdir='wspd')
   # clwp_values = clwp.calculate_lwp(filein=files)
   # uv.values = clwp_values
    uv_xr = wrf_assign_coords(uv.metpy.assign_crs(
        grid_mapping_name="latitude_longitude", earth_radius=6371229.0
    ))
    uv_xr = uv_xr.rename({"south_north": "y", "west_east": "x"})
    uv_xr.Time.values = uv.Time.values
    uv_cs_pst.append(cross_section(uv_xr, start, end))




plt.figure(figsize=(12, 6))
(xr.concat(uv_cs_pst, dim='Time').T-xr.concat(uv_cs_pre, dim='Time').T).plot(cmap='coolwarm')
plt.show()





