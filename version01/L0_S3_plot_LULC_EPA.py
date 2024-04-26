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


plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


domain='3'
geog_file_pattern = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/update_geog/def_geog_files/geo_em.d0{domain}_epa_{year}.nc'

results = pd.DataFrame(columns=[f'LULC_Class_{i}' for i in range(1, 18)])

for year in range(2030, 2101, 10):
    geog_file = geog_file_pattern.format(year=year, domain=domain)
    
    pre_geog_file = xr.open_dataset(geog_file, engine="netcdf4")
    wrf_lulc_pre = pre_geog_file["LU_INDEX"].squeeze().squeeze().values
    
    wrf_lulc_pre_count = np.array([
        np.where(wrf_lulc_pre == lulc_classes_id)[0].shape[0]
        for lulc_classes_id in np.arange(1, 18)
    ])
    
    results.loc[f'{year}'] = wrf_lulc_pre_count

results.columns=list(get_lulc_colormap()[1].keys())
results = results.T

results.to_csv(
    "../tables/EPA_" + ".".join(geog_file.split("/")[-1].split(".")[:2]) + "_LULC_count.csv"
)

plt.figure(figsize=(12, 5))
plt.plot(list(results.loc['Urban and Built-Up'].index), results.loc['Urban and Built-Up'].values, 'r-o')
plt.xlabel('Years')
plt.ylabel(f'Urban Pixel Count (d0{domain})')
plt.tight_layout()
plt.savefig("../tables/EPA_" + ".".join(geog_file.split("/")[-1].split(".")[:2]) + "_LULC_count.jpeg", dpi=400)
#plt.show()


