import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from src.coast import plot_coast, plot_houston
import cartopy.crs as ccrs

plt.rcParams.update({"font.size": 14, "font.weight": "bold", "savefig.dpi": 300})

geog_file_old = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/pre/WPS/geo_em.d01.nc"
geog_file_new = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post/WPS/geo_em.d01.nc"



A = xr.open_dataset(geog_file_old, engine="netcdf4")
wrf_longitudes = A["XLONG_M"].squeeze().values
wrf_latitudes = A["XLAT_M"].squeeze().values
wrf_lulc_old = A["LU_INDEX"].squeeze().squeeze().values
wrf_grnf_old = A["GREENFRAC"].sel(month=10).squeeze().values

A = xr.open_dataset(geog_file_new, engine="netcdf4")
wrf_longitudes = A["XLONG_M"].squeeze().values
wrf_latitudes = A["XLAT_M"].squeeze().values
wrf_lulc_new = A["LU_INDEX"].squeeze().squeeze().values
wrf_grnf_new = A["GREENFRAC"].sel(month=10).squeeze().values


fig, axs = plt.subplots(1, 2, figsize=(16, 6), subplot_kw={'projection': ccrs.PlateCarree()})
out1 = axs[0].contourf(
    wrf_longitudes, wrf_latitudes, wrf_grnf_old, levels=np.arange(0, 1, 0.1), cmap="BuGn", extend='both'
)

out2 = axs[1].contourf(
    wrf_longitudes, wrf_latitudes, wrf_grnf_new, levels=np.arange(0, 1, 0.1), cmap="BuGn", extend='both'
)


[plot_coast(axs) for axs in axs]
plt.colorbar(out1, ax=axs[0], shrink=0.8)
plt.colorbar(out2, ax=axs[1], shrink=0.8)
plot_houston(axs[0])
plot_houston(axs[1])

axs[0].axes.set_aspect("equal")
axs[1].axes.set_aspect("equal")

axs[0].set_title('GREENFRAC (LULC 2001)')
axs[1].set_title('GREENFRAC (LULC 2017)')
plt.tight_layout()
#plt.savefig('.'.join(geog_file_old.split("/")[-1].split(".")[:2])+'_GREENFRAC_map.jpeg', dpi=300)

plt.show()



