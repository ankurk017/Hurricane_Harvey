from src.wrf_src import wrf_assign_coords, find_common_min_max
import progressbar
import xarray as xr
from wrf import getvar
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
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
        if mean:
            crop_values = var_pre.sel(west_east=slice(location[0]-box, location[0]+box), south_north=slice(location[1]-box, location[1]+box)).mean()
        else:
            crop_values = var_pre.sel(west_east=slice(location[0]-box, location[0]+box), south_north=slice(location[1]-box, location[1]+box))

        var_timeseries_pre_f.append(crop_values)
    var_descrip = f"{var_pre.attrs['description']} ({var_pre.attrs['units']})"
    return xr.concat(var_timeseries_pre_f, dim="Time"), var_descrip


home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/PBL/'

wrfoutfile_pbl00 = sorted(glob.glob(home_2512 + f'/WRF_PBL00//test/em_real/wrfout_d02_2017-*'))
wrfoutfile_pbl01 = sorted(glob.glob(home_2512 + f'/WRF_PBL01//test/em_real/wrfout_d02_2017-*'))
wrfoutfile_pbl02 = sorted(glob.glob(home_2512 + f'/WRF_PBL02//test/em_real/wrfout_d02_2017-*'))
wrfoutfile_pbl04 = sorted(glob.glob(home_2512 + f'/WRF_PBL04//test/em_real/wrfout_d02_2017-*'))

location = (-97.061, 27.8339)  # Landfall location
location = (-95.499,  29.74)  # Houston location
box=0.75

variable_name = 'PBLH'
variable_name = 'UST'
variable_name = 'ACHFX'
variable_name = 'ACLHF'


var_timeseries_pbl00, var_descrip = get_wrf_var(wrfoutfile_pbl00, location, variable_name = variable_name, box=box)
var_timeseries_pbl01, _ = get_wrf_var(wrfoutfile_pbl01, location, variable_name = variable_name, box=box)
var_timeseries_pbl02, _ = get_wrf_var(wrfoutfile_pbl02, location, variable_name = variable_name, box=box)
var_timeseries_pbl04, _ = get_wrf_var(wrfoutfile_pbl04, location, variable_name = variable_name, box=box)

labels = ['No PBL', 'YSU', 'MYJ TKE', 'EDMF QNSE']
fig, axs = plt.subplots(1, 1, figsize=(10, 5.5))

axs.plot(var_timeseries_pbl04["Time"],  var_timeseries_pbl00.mean(dim=['west_east', 'south_north']),    "b-",    label=labels[0])
axs.plot(var_timeseries_pbl04["Time"],  var_timeseries_pbl01.mean(dim=['west_east', 'south_north']),    "r-",    label=labels[1])
axs.plot(var_timeseries_pbl04["Time"],  var_timeseries_pbl02.mean(dim=['west_east', 'south_north']),    "g-",    label=labels[2])
axs.plot(var_timeseries_pbl04["Time"],  var_timeseries_pbl04.mean(dim=['west_east', 'south_north']),    "m-",    label=labels[3])

axs.set_xlabel("Time")
axs.set_ylabel(var_descrip)
plt.legend()
plt.xticks(rotation=25)
plt.xlim(
    (var_timeseries_pbl04["Time"][0], var_timeseries_pbl04["Time"][-1])
)

plt.grid()
plt.minorticks_on()
axs.xaxis.set_minor_locator(ticker.MaxNLocator(46))
plt.tight_layout()
plt.savefig(f'../figures_paper/pbl/timeseries_{variable_name}.jpeg')


data_matrices = [val.sel(Time='2017-08-27').mean(dim='Time') for val in [var_timeseries_pbl01, var_timeseries_pbl02, var_timeseries_pbl04]]
min_val, max_val = find_common_min_max(data_matrices)

fig, axes = plt.subplots(1, 3, figsize=(15, 6), subplot_kw={'projection': ccrs.PlateCarree()}, sharey=True, gridspec_kw={'bottom': 0.2, 'left':0.09, 'right':0.96, 'wspace':0.3})

for ax, data_matrix, label in zip(axes.ravel(), data_matrices, labels[1:]):
    pout = ax.contourf(data_matrix.west_east, data_matrix.south_north, data_matrix, levels=np.linspace(min_val, max_val, 50), extend='both', cmap='jet')
    coast.plot_coast(ax, houston=True, houston_color='black', houston_linewidth=0.6)
    ax.set_title(label)

cax = fig.add_axes([0.15, 0.12, 0.7, 0.042])  
cbar = plt.colorbar(pout, cax=cax, orientation='horizontal', shrink=0.6)  

cbar.set_label(var_descrip, labelpad=10)  
plt.tight_layout()
plt.savefig(f'../figures_paper/pbl/spatial_{variable_name}.jpeg')
#plt.show()




