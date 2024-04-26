from src.wrf_src import wrf_assign_coords, find_common_min_max
import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords, interplevel
#from self_utils import ahps, coast
import src.coast as coast
import numpy as np
import cartopy.io.shapereader as shpreader
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.prepared import prep
import datetime
from matplotlib.patches import Rectangle
from matplotlib import ticker 


plt.rcParams.update({"font.size": 16, "font.weight": "bold"})

#		UST:description = "U* IN SIMILARITY THEORY" ;
#		PBLH:description = "PBL HEIGHT" ;

def get_wrf_var(wrfoutfile_pre_f, location_f, variable_name = 'UST', box=0.75, mean=False):
    var_timeseries_pre_f = []

    for wrf_files_id in progressbar.progressbar(range(len(wrfoutfile_pre_f) - 1)):
        wrf_ncfile_pre = Dataset(wrfoutfile_pre_f[wrf_files_id + 1])

        var_pre = wrf_assign_coords(getvar(wrf_ncfile_pre, variable_name))

        z = getvar(wrf_ncfile_pre, "z", units="km")
        wspd = getvar(wrf_ncfile_pre, "wspd_wdir", units="m/s")[0,:]
        wspd_1 = interplevel(wspd, z, 0.75)
        wspd_2 = interplevel(wspd, z, 3.25)
        var_pre = wrf_assign_coords(wspd_1-wspd_2)

        if mean:
            crop_values = var_pre.sel(west_east=slice(location[0]-box, location[0]+box), south_north=slice(location[1]-box, location[1]+box)).mean()
        else:
            crop_values = var_pre.sel(west_east=slice(location[0]-box, location[0]+box), south_north=slice(location[1]-box, location[1]+box))

        var_timeseries_pre_f.append(crop_values)
    #var_descrip = f"{var_pre.attrs['description']} ({var_pre.attrs['units']})"
    var_descrip = "Wind Shear (m/s)"
    return xr.concat(var_timeseries_pre_f, dim="Time"), var_descrip



home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
wrfoutfile_pbl02 = sorted(
    glob.glob(home_2512 + f"/pre/WRF/test/em_real/wrfout_d02_2017-08-27*") + glob.glob(home_2512 + f"/pre/WRF/test/em_real/wrfout_d02_2017-08-28*")
)[:24]
wrfoutfile_pbl04 = sorted(
    glob.glob(home_2512 + f"/post/WRF/test/em_real/wrfout_d02_2017-08-27*") + glob.glob(home_2512 + f"/post/WRF/test/em_real/wrfout_d02_2017-08-28*")
)[:24]


location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location

box=0.15


locations = [(-95., 30.01), (-94.6, 29.92), (-95.71, 28.8), (-96.72, 28.90)]
labels = ['R1 @ 07 UTC', 'R1 @ 12 UTC', 'R2 @ 12 UTC (coast)', 'R3 @ 08 UTC']
variable_name = 'T2' # 'PBLH'  # or 'T2'

fig, axs = plt.subplots(2, 2, figsize=(15, 8), sharex=True, sharey=True)

for i, location in enumerate(locations):
    var_timeseries_pbl02, var_descrip = get_wrf_var(wrfoutfile_pbl02, location, variable_name=variable_name, box=box)
    var_timeseries_pbl04, _ = get_wrf_var(wrfoutfile_pbl04, location, variable_name=variable_name, box=box)

    row = i // 2
    col = i % 2

    axs[row, col].plot(var_timeseries_pbl04["Time"], var_timeseries_pbl02.mean(dim=['west_east', 'south_north']), "b-o", label='LULC 2001')
    axs[row, col].plot(var_timeseries_pbl04["Time"], var_timeseries_pbl04.mean(dim=['west_east', 'south_north']), "r-o", label='LULC 2017')

    axs[row, col].set_title(labels[i])
    axs[row, col].grid()
    axs[row, col].minorticks_on()
    axs[row, col].xaxis.set_minor_locator(ticker.MaxNLocator(46))

axs[1, 0].set_xlabel("Time")
axs[1, 1].set_xlabel("Time")
axs[0, 0].set_ylabel(var_descrip)
axs[1, 0].set_ylabel(var_descrip)

for ax in axs.flat:
    ax.xaxis.set_tick_params(rotation=25)

plt.legend()
plt.xlim((var_timeseries_pbl04["Time"][0], var_timeseries_pbl04["Time"][-1]))

plt.tight_layout()
#plt.savefig(f'../figures_paper/cold_pools_timeseries_radius_{int(box*100)}deg.jpeg')
plt.show()



