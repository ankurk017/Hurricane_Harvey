import tropycal.tracks as tracks
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from wrf import getvar, interplevel, latlon_coords
import cartopy.crs as ccrs

import src.coast

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


def get_wrfvar(wrf_ncfile, variable_name, pres_level):
    pressure_3d = getvar(wrf_ncfile, "pressure")
    variable = xr.concat(
        [interplevel(getvar(wrf_ncfile, variable_name),
                     pressure_3d, pres_level)],
        dim="level",
    )
    return variable


def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"


basin = tracks.TrackDataset(basin="north_atlantic",
                            source="hurdat", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))


data_dir = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_IDA/Analysis/wrfouts/"
data_dir = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/"
data_dir = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"


pres_levels = np.array(
    (np.arange(1000, 700, -25), np.arange(700, 100, -50))).ravel()

wrf_reffile = data_dir + "wrfout_d02_2017-08-24_06:00:00"
ref_ncfile = Dataset(wrf_reffile)

ref_theta = get_wrfvar(ref_ncfile, "theta", pres_levels)


wrf_file = data_dir + "wrfout_d02_2017-08-26_00:00:00"
wrf_file = data_dir + "wrfout_d02_2017-08-27_12:00:00"

ncfile = Dataset(wrf_file)

slp = getvar(ncfile, "slp")
lats, lons = latlon_coords(slp)
xaxis = (
    lons.mean(dim="south_north") - lons.mean(dim="south_north").mean()
).values * 111.11


theta = get_wrfvar(ncfile, "theta", pres_levels)
rh = get_wrfvar(ncfile, "rh", pres_levels)


fig, ax = plt.subplots(
    1,
    2,
    figsize=(13, 5.5),
    subplot_kw={"projection": ccrs.PlateCarree()},
    sharex=True,
    sharey=True,
)

u_wind = (
    get_wrfvar(ncfile, "ua", [850])
    .isel(
        south_north=np.arange(
            129 - 70,
            129 + 70,
        ),
        west_east=np.arange(
            129 - 70,
            129 + 70,
        ),
    )
    .squeeze()
)
v_wind = (
    get_wrfvar(ncfile, "va", [850])
    .isel(
        south_north=np.arange(
            129 - 70,
            129 + 70,
        ),
        west_east=np.arange(
            129 - 70,
            129 + 70,
        ),
    )
    .squeeze()
)
wind_upper = np.round(u_wind.max().values/10)*20
wind_upper = 80

ax[0].contourf(
    u_wind["XLONG"],
    u_wind["XLAT"],
    np.sqrt(u_wind**2 + v_wind**2),
    cmap="bwr",
    levels=np.arange(0, wind_upper, 5),
)
ax[0].quiver(
    u_wind["XLONG"][::10, ::10],
    u_wind["XLAT"][::10, ::10],
    u_wind[::10, ::10],
    v_wind[::10, ::10],
)

src.coast.plot_coast(ax[0])
ax[0].plot(harvey["lon"], harvey["lat"], "b-")
ax[0].set_xlim((u_wind["XLONG"].min().values, u_wind["XLONG"].max().values))
ax[0].set_ylim((u_wind["XLAT"].min().values, u_wind["XLAT"].max().values))
ax[0].set_title(f"{str(slp['Time'].values)[:13]} | 850 hPa")
ax[0].plot(getvar(ncfile, 'slp')['XLONG'][np.where(getvar(ncfile, 'slp') == getvar(ncfile, 'slp').min())].values,
           getvar(ncfile, 'slp')['XLAT'][np.where(getvar(ncfile, 'slp') == getvar(ncfile, 'slp').min())].values, 'ro')

u_wind = (
    get_wrfvar(ncfile, "ua", [350])
    .isel(
        south_north=np.arange(
            129 - 70,
            129 + 70,
        ),
        west_east=np.arange(
            129 - 70,
            129 + 70,
        ),
    )
    .squeeze()
)
v_wind = (
    get_wrfvar(ncfile, "va", [350])
    .isel(
        south_north=np.arange(
            129 - 70,
            129 + 70,
        ),
        west_east=np.arange(
            129 - 70,
            129 + 70,
        ),
    )
    .squeeze()
)

cont = ax[1].contourf(
    u_wind["XLONG"],
    u_wind["XLAT"],
    np.sqrt(u_wind**2 + v_wind**2),
    cmap="bwr",
    levels=np.arange(0, wind_upper, 5),
)

ax[1].quiver(
    u_wind["XLONG"][::10, ::10],
    u_wind["XLAT"][::10, ::10],
    u_wind[::10, ::10],
    v_wind[::10, ::10],
)
src.coast.plot_coast(ax[1])
ax[1].plot(harvey["lon"], harvey["lat"], "b-")
ax[1].set_xlim((u_wind["XLONG"].min().values, u_wind["XLONG"].max().values))
ax[1].set_ylim((u_wind["XLAT"].min().values, u_wind["XLAT"].max().values))
ax[1].set_title(f"{str(slp['Time'].values)[:13]} | 350 hPa")
ax[1].plot(getvar(ncfile, 'slp')['XLONG'][np.where(getvar(ncfile, 'slp') == getvar(ncfile, 'slp').min())].values,
           getvar(ncfile, 'slp')['XLAT'][np.where(getvar(ncfile, 'slp') == getvar(ncfile, 'slp').min())].values, 'ro')

plt.tight_layout()
cbar = plt.colorbar(cont, ax=ax.ravel())
cbar.set_label("Wind Speed (m/s)")
plt.savefig(f"../figures/Quiver_ws_{str(slp['Time'].values)[:13]}.jpeg")
# plt.show()

theta_anom = theta - ref_theta

theta_anom_data = theta_anom.isel(
    south_north=np.arange(
        129 - 70,
        129 + 70,
    )
).mean(dim="south_north")

rh_data = rh.isel(
    south_north=np.arange(
        129 - 70,
        129 + 70,
    )
).mean(dim="south_north")


plt.figure(figsize=(6.5, 5.5))
theta_plot = plt.contourf(
    xaxis,
    pres_levels,
    theta_anom_data,
    levels=np.arange(-10, 11, 1),
    cmap="coolwarm",
)
rh_plot = plt.contour(
    xaxis,
    pres_levels,
    rh_data,
    6,
    colors="k",
)
plt.clabel(rh_plot, fontsize=14, inline=1, fmt=fmt)

tilt = np.where(
    theta_anom_data
    == theta_anom_data.isel(west_east=np.arange(129 - 70, 129 + 70)).max(
        dim="west_east"
    )
)
plt.plot(xaxis[tilt[1]], pres_levels[tilt[0]], "m--")
plt.plot([0, 0], [pres_levels.min(), pres_levels.max()], 'k--')
plt.gca().invert_yaxis()
cbar = plt.colorbar(theta_plot)
cbar.set_label("Potential Temperature Anomaly (K)")
plt.xlabel("Radius of the storm center (km)")
plt.ylabel("Pressure (hPa)")
plt.title(f"{str(slp['Time'].values)[:13]} ")
plt.tight_layout()
plt.savefig(f"../figures/Theta_RH_{str(slp['Time'].values)[:13]}.jpeg")
plt.show()
