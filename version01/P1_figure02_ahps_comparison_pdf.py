from sklearn.metrics import mean_squared_error
import cartopy.io.shapereader as shpreader
from matplotlib import ticker, cm
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
from netCDF4 import Dataset
from wrf import getvar, latlon_coords

# from self_utils import ahps, coast
import src.coast as coast
import src.ahps as ahps
import numpy as np
from scipy.interpolate import griddata
from scipy import stats


plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


def rmse(predictions, targets):
    return np.sqrt(np.nanmean((predictions - targets) ** 2))


def correlation_without_nans(array1, array2):
    # Remove NaN values from both arrays
    valid_indices = np.logical_and(~np.isnan(array1), ~np.isnan(array2))
    valid_array1 = array1[valid_indices]
    valid_array2 = array2[valid_indices]

    # Calculate correlation coefficient
    correlation = np.corrcoef(valid_array1, valid_array2)[0, 1]

    return correlation


def calc_wrf_pcp_error(date, case, domain_bb=[-100, -93, 25.5, 31.5]):
    home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/AHPS/"
    ahps_files = glob.glob(home + "nws_precip*_conus.nc")
    home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/"

    home_2512 = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/"
    wrfoutfile = sorted(
        glob.glob(home_2512 + f"/{case}/WRF_cntl/test/em_real/wrfout_d02_2017-*{date}*")
    )

    wrf_pcp = (
        getvar(Dataset(wrfoutfile[-1]), "RAINC")
        + getvar(Dataset(wrfoutfile[-1]), "RAINNC")
    ) - (
        getvar(Dataset(wrfoutfile[0]), "RAINC")
        + getvar(Dataset(wrfoutfile[0]), "RAINNC")
    )
    slp = getvar(Dataset(wrfoutfile[0]), "slp")
    wrf_lat, wrf_lon = latlon_coords(slp)

    ahps_fileid = np.where([f"201708{date}" in files for files in ahps_files])[0][0]
    ahps_lon, ahps_lat, ahps_pcp = ahps.read_ahps(ahps_files[ahps_fileid])

    merged_lats = np.arange(domain_bb[2], domain_bb[3], 0.05)
    merged_lons = np.arange(domain_bb[0], domain_bb[1], 0.05)
    merged_meshlon, merged_meshlat = np.meshgrid(merged_lons, merged_lats)

    ahps_pcp_domain = griddata(
        (ahps_lon.ravel(), ahps_lat.ravel()),
        ahps_pcp.ravel(),
        (merged_meshlon, merged_meshlat),
    )
    wrf_pcp_domain = griddata(
        (wrf_lon.values.ravel(), wrf_lat.values.ravel()),
        wrf_pcp.values.ravel(),
        (merged_meshlon, merged_meshlat),
    )
    wrf_pcp_error = wrf_pcp_domain.ravel() - ahps_pcp_domain.ravel()
    # rms = rmse(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
    # cor = correlation_without_nans(ahps_pcp_domain.ravel(), wrf_pcp_domain.ravel())
    return wrf_pcp_error


########## plotting differenbce

from scipy.stats import gaussian_kde

fig, ax = plt.subplots(2, 2, figsize=(14, 9))
plt.subplots_adjust(hspace=0.5)


for ddate, iindex in zip(("27", "28", "29", "30"), (0, 1, 2, 3)):
    print(ddate)
    dom = [-97, -94.5, 28, 31]
    pre_error = calc_wrf_pcp_error(ddate, "pre", domain_bb=dom)
    post_error = calc_wrf_pcp_error(ddate, "post", domain_bb=dom)

    pdata = (pre_error, post_error)
    col = ("b", "r")
    lab = ("LULC 2001", "LULC 2017")


    for index in range(len(pdata)):
        data = pdata[index][~np.isnan(pdata[index])]

        kde = gaussian_kde(data)
        x = np.linspace(np.min(data), np.max(data), 100)
        pdf = kde(x)
        ax.ravel()[iindex].plot(x, pdf, label=lab[index], color=col[index], linewidth=3)
        ax.ravel()[iindex].axvline(
            x=x[pdf == pdf.max()], color=col[index], linestyle="--", linewidth=2
        )

    ax.ravel()[iindex].axvline(x=0, color="black", linestyle="--")

    ax.ravel()[iindex].set_title(f"201708{ddate}")
    ax.ravel()[iindex].set_xlabel(" WRF precipitation error (mm/hr)")
    ax.ravel()[iindex].set_ylabel("Probability")

    ax.ravel()[iindex].legend()
    plt.tight_layout()
plt.savefig("../figures_paper/ahps_pdf/pdf_ahps.jpeg", dpi=300)
plt.show()

