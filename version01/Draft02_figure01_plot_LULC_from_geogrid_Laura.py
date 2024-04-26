import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import matplotlib.patches as patches
from osgeo import gdal, osr
import pyproj

from scipy.interpolate import griddata
from src.lulc_colormap import get_lulc_colormap
from src.coast import plot_coast


plt.rcParams.update({"font.size": 18, "font.weight": "bold"})


def get_bb(file_name):
    return np.array(
        (
            file_name["XLONG_M"].min().values,
            file_name["XLONG_M"].max().values,
            file_name["XLAT_M"].min().values,
            file_name["XLAT_M"].max().values,
        )
    )


def draw_box(bb):
        return patches.Rectangle(
            (bb[0], bb[2]),
            bb[1] - bb[0],
            bb[3] - bb[2],
            linewidth=4,
            edgecolor="m",
            facecolor="none",
            alpha=1,
        )



def plot_lulc_geogrid(geog_file, label='LULC 2001', legend=False, axes=None):
    pre_geog_file = xr.open_dataset(geog_file, engine="netcdf4")

    wrf_longitudes_pre = pre_geog_file["XLONG_M"].squeeze().values
    wrf_latitudes_pre = pre_geog_file["XLAT_M"].squeeze().values
    wrf_lulc_pre = pre_geog_file["LU_INDEX"].squeeze().squeeze().values

    lulc_cmap, lulc_classes = get_lulc_colormap()


    img = axes.pcolormesh(
        wrf_longitudes_pre,
        wrf_latitudes_pre,
        wrf_lulc_pre,
        vmin=1,
        vmax=18,
        cmap=lulc_cmap,
        shading="auto",
    )
    plot_coast(axes)
    if legend:
     cbar = plt.colorbar(img, ticks=np.arange(1, 18, 1) + 0.5, shrink=0.8)
     cbar.ax.set_yticklabels(list(lulc_classes.keys()), fontsize=12)  # horizontal colorbar
     axes.set_extent(get_bb(pre_geog_file))
    urban_count = f'{np.where(wrf_lulc_pre==13)[0].shape[0]} ({np.round(np.where(wrf_lulc_pre==13)[0].shape[0]*100/(wrf_lulc_pre.ravel().shape[0]), 2)} %)'

    axes.set_title(f'{label} | {urban_count}')
    axes.set_xlim((-96.60, -92.47))
    axes.set_ylim((28.55, 31.17))
    #plot_domain(axes)
    plt.tight_layout()
    #plt.savefig('.'.join(geog_file.split("/")[-1].split(".")[:2])+'_LULC_pre_map.jpeg', dpi=300)
    return axes

def plot_domain(ax):
 geog_files = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WPS/geo_em.d*.nc'))
 d1bb = get_bb(xr.open_dataset(geog_files[0]))
 d2bb = get_bb(xr.open_dataset(geog_files[1]))
 d3bb = get_bb(xr.open_dataset(geog_files[2]))
 [ax.add_patch(draw_box(dbb)) for dbb in (d1bb, d2bb, d3bb)]
 return None

domain = 'd02'
geog_2001 = f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WPS/geo_em.{domain}.nc"
geog_2020 = f"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WPS/geo_em.{domain}.nc"
geog_2050 = f'/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/update_geog/def_geog_files/geo_em.{domain}_epa_2050_cropland.nc'
geog_2100 = f'/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/update_geog/def_geog_files/geo_em.{domain}_epa_2100_cropland.nc'

fig = plt.figure(figsize=(12.5, 9))

# Create subplots with specific projection
axes1 = plt.subplot(2, 2, 1, projection=ccrs.PlateCarree())
axes2 = plt.subplot(2, 2, 2, projection=ccrs.PlateCarree())
axes3 = plt.subplot(2, 2, 3, projection=ccrs.PlateCarree())
axes4 = plt.subplot(2, 2, 4, projection=ccrs.PlateCarree())

plot_lulc_geogrid(geog_2001, label='(a) LULC 2001', axes=axes1)
plot_lulc_geogrid(geog_2020, label='(b) LULC 2020', axes=axes2)
plot_lulc_geogrid(geog_2050, label='(c) LULC 2050', axes=axes3)
plot_lulc_geogrid(geog_2100, label='(d) LULC 2100', axes=axes4)

plt.subplots_adjust(hspace=0.25)
#plt.savefig(f'../figures_draft01/Laura_fig01_{domain}_cropland.jpeg', dpi=400)
plt.show()

