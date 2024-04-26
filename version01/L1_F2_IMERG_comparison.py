from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
from self_utils import ahps, coast
import numpy as np


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_IDA/Analysis/"

ahps_files = glob.glob(home + "AHPS/nws_precip*_conus.nc")


wrfoutfile = sorted(glob.glob("../../wrfouts/wrfout_d01*-08-30*_*:00:00"))
wrf_pcp = (
    getvar(Dataset(wrfoutfile[-1]), "RAINC") +
    getvar(Dataset(wrfoutfile[-1]), "RAINNC")
) - (
    getvar(Dataset(wrfoutfile[0]), "RAINC") +
    getvar(Dataset(wrfoutfile[0]), "RAINNC")
)
slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)

ahps_fileid = np.where(['20210830' in files for files in ahps_files])[0][0]
ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])


domain_bb = [-95, -85, 27.5, 35]

fig, ax = plt.subplots(
    1, 1, figsize=(11, 6), subplot_kw={"projection": ccrs.PlateCarree()}
)
cont = ax.contourf(ahps_lon, ahps_lat, ahps_pcp, cmap="jet")
plt.colorbar(cont)
coast.plot_coast(ax)
plt.xlim((domain_bb[0], domain_bb[1]))
plt.ylim((domain_bb[2], domain_bb[3]))

fig, ax = plt.subplots(
    1, 1, figsize=(11, 6), subplot_kw={"projection": ccrs.PlateCarree()}
)
cont = ax.contourf(wrf_lon, wrf_lat, wrf_pcp, cmap="jet")
plt.colorbar(cont)
coast.plot_coast(ax)
plt.xlim((domain_bb[0], domain_bb[1]))
plt.ylim((domain_bb[2], domain_bb[3]))
plt.show()

merged_lats = np.arange(domain_bb[2], domain_bb[3], 0.05)
merged_lons = np.arange(domain_bb[0], domain_bb[1], 0.05)
merged_meshlon, merged_meshlat = np.meshgrid(merged_lons, merged_lats)


ahps_pcp_domain = griddata((ahps_lon.ravel(), ahps_lat.ravel(
)), ahps_pcp.ravel(), (merged_meshlon, merged_meshlat))
wrf_pcp_domain = griddata((wrf_lon.values.ravel(), wrf_lat.values.ravel(
)), wrf_pcp.values.ravel(), (merged_meshlon, merged_meshlat))
