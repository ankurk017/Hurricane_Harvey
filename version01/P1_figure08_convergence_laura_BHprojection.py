import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
#from self_utils import ahps, coast
import src.coast as coast
import numpy as np
from scipy.interpolate import griddata
import progressbar
import matplotlib
from src.wrf_src import wrf_assign_coords, crop_region
from metpy.units import units
import metpy.calc as mpcalc
from wrf import interplevel

# matplotlib.use('Agg')

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


domain_bb = [-100, -93, 25.5, 31.5]
domain_bb = [-96.72, -94.21, 28.71, 30.91]
urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 1.2}

urban_change = {"south_north": 29.74, "west_east": -95.49, "box": 2.65}
urban_change = {"south_north": 30.00, "west_east": -94.49, "box": 2.65}

domain_bb = [urban_change['west_east']-urban_change['box'], urban_change['west_east']+urban_change['box'], urban_change['south_north']-urban_change['box'], urban_change['south_north']+urban_change['box']]


levels  = np.linspace(-1e-3, 1e-3, 10)
levels  = np.linspace(-1, 1, 12)
diff_levels = np.arange(-50, 55, 5)

home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Laura1/"

wrfoutfile_2001 = sorted(glob.glob('/rgroup/stela/akumar/BH_Simulations/2000/WRF_452_BH/test/em_real/wrfout_d02*'))[13:43]
wrfoutfile_2017 = sorted(glob.glob('/rgroup/stela/akumar/BH_Simulations/2020/WRF_452_BH/test/em_real/wrfout_d02*'))[13:43]
wrfoutfile_2050 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD//WRF_Building_Heights/WRF_Simulations/WRF_2050/WRF/test/em_real/wrfout_d02*'))[13:43]
wrfoutfile_2100 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD//WRF_Building_Heights/WRF_Simulations/WRF_2100/WRF/test/em_real/wrfout_d02*'))[13:43]


prefiles = wrfoutfile_2001
postfiles = wrfoutfile_2017

years = (2001, 2017, 2050, 2100)


for fileid in progressbar.progressbar(range(8, len(prefiles) - 1)):

    wrf_pcp_pre_pres = getvar(Dataset(prefiles[fileid + 1]), "pres", units='hPa')
    wrf_pcp_pre_uv = crop_region(wrf_assign_coords(interplevel(getvar(Dataset(prefiles[fileid + 1]), "uvmet"), wrf_pcp_pre_pres, 750)), location=urban_change)

    wrf_pcp_post_pres = getvar(Dataset(postfiles[fileid + 1]), "pres", units='hPa')
    wrf_pcp_post_uv = crop_region(wrf_assign_coords(interplevel(getvar(Dataset(postfiles[fileid + 1]), "uvmet"), wrf_pcp_post_pres, 750)), location=urban_change)

    divergence_pre = mpcalc.divergence(wrf_pcp_pre_uv.sel(u_v="u").values, wrf_pcp_pre_uv.sel(u_v="v").values, dx=units.Quantity(0.027*111.11*1000, "m"), dy=units.Quantity(0.027*111.11*1000, "m"))*1e3
    divergence_post = mpcalc.divergence(wrf_pcp_post_uv.sel(u_v="u").values, wrf_pcp_post_uv.sel(u_v="v").values, dx=units.Quantity(0.027*111.11*1000, "m"), dy=units.Quantity(0.027*111.11*1000, "m"))*1e3

    slp = getvar(Dataset(prefiles[fileid + 1]), "slp")
    wrf_lat, wrf_lon = latlon_coords(slp)
#    pre_post = (wrf_pcp_pre, wrf_pcp_post)
    pre_post_uv = (wrf_pcp_pre_uv, wrf_pcp_post_uv)
    pre_post_div = (divergence_pre, divergence_post)

    levels  = np.linspace(-0.25, 0.25, 9)
    fig, ax = plt.subplots(   2,    1,     figsize=(6.7, 11),        sharey=True,  sharex=True,       subplot_kw={"projection": ccrs.PlateCarree()}, 
 gridspec_kw={'wspace': 0.4, })


    for pre_post_id in range(2):
        cont = ax[pre_post_id].contourf(
            pre_post_uv[pre_post_id]["west_east"],
            pre_post_uv[pre_post_id]["south_north"],
            pre_post_div[pre_post_id],
            cmap="PiYG",
            levels=levels, extend='both', 
        )


        ax[pre_post_id].quiver(
            pre_post_uv[pre_post_id]["west_east"][::5],
            pre_post_uv[pre_post_id]["south_north"][::5],
            pre_post_uv[pre_post_id].sel(u_v="u")[::5, ::5] ,
            pre_post_uv[pre_post_id].sel(u_v="v")[::5, ::5],  units='x', scale=110, width=0.019, 
        )


        ax[pre_post_id].set_xlim((domain_bb[0], domain_bb[1]))
        ax[pre_post_id].set_ylim((domain_bb[2], domain_bb[3]))
        ax[pre_post_id].set_title(f"WRF (LULC {years[pre_post_id]})")
        ax[0].set_xticklabels([])

        ax[pre_post_id].tick_params(axis='x', labelrotation=45)

        if pre_post_id==1:
            cbar = plt.colorbar(cont, ax=ax[:2].ravel(), ticks=levels, shrink=0.75)
            cbar.ax.set_yticklabels(np.round(levels, 3))
            cbar.ax.set_ylabel(r"Divergence x $10^{-3}$ (1/s)")
       	coast.plot_coast(ax[pre_post_id], gridlines_alpha=0.1, houston=True, houston_linewidth=0.7)
    """
    cont = ax[2].contourf(pre_post_uv[pre_post_id]["west_east"], pre_post_uv[pre_post_id]["south_north"], divergence_post-divergence_pre, cmap="bwr")
    coast.plot_coast(ax[2], gridlines_alpha=0, houston=True)
    ax[2].set_xlim((domain_bb[0], domain_bb[1]))
    ax[2].set_ylim((domain_bb[2], domain_bb[3]))
    ax[2].set_title(f'WRF (LULC {years[pre_post_id]})')
    cbar = plt.colorbar(cont, ax=ax[2], shrink=0.75)
    # 	cbar = plt.colorbar(cont, ax=ax[2], ticks=diff_levels, shrink=0.75)
    # 	cbar.ax.set_yticklabels(diff_levels)
    # 	cbar.ax.set_ylabel('Precipitation Error (mm)')
#    plt.tight_layout()
    """
    plt.suptitle(str(slp.Time.values)[:13], fontweight="bold")
    plt.savefig(f"../figures/Laura/Ensemble/Divergence//Houston_divergence_{str(slp.Time.values)[:13]}.jpeg", dpi=400)
    plt.close()

#plt.show()





