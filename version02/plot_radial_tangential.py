import matplotlib.pyplot as plt
import glob
import xarray as xr
import numpy as np

plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "font.family": "monospace"}
)


def plot_rad_tang(uv_cs, xx="radius", case_name='pre', prefix=None, row = 0):

    # Radial Plot
    rad = axs[row, 0].contourf(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        cmap="bwr",
        levels=np.arange(-30, 35, 5),
    )
    contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(-60, 60, 10),
    )
    neg_contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(-20, 0, 5),
    )
    pos_contour_rad = axs[row, 0].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["radial"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, 40, 10),
    )
    axs[row, 0].clabel(contour_rad, inline=True, fmt="%1.1f", fontsize=12)

    skip = {'x': 4, 'y':1}

    rad_vals = uv_cs['wa_interp']
    xxx, yyy, = np.meshgrid(uv_cs[xx][::skip['x']], uv_cs['level'][::skip['y']])
    axs[row, 0].scatter(xxx, yyy, rad_vals[::skip['y'], ::skip['x']]*100, c=rad_vals[::skip['y'], ::skip['x']]*70, alpha=0.50, cmap='Greens', marker='o')
    axs[row, 0].scatter(xxx, yyy, -(rad_vals[::skip['y'], ::skip['x']]*500), c=rad_vals[::skip['y'], ::skip['x']]*70, alpha=1, cmap='Reds', marker='x')

    # Tangential Wind
    tan = axs[row, 1].contourf(uv_cs[xx], uv_cs["level"], uv_cs["tangential"], cmap="jet", levels=np.arange(-15, 55, 5))
    contour_tan = axs[row, 1].contour(
        uv_cs[xx], uv_cs["level"], uv_cs["tangential"], colors="k", linewidths=0.5
    )
    neg_contour_tan = axs[row, 1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        linestyles="dashed",
        levels=np.arange(-20, 0, 5),
    )
    pos_contour_tan = axs[row, 1].contour(
        uv_cs[xx],
        uv_cs["level"],
        uv_cs["tangential"],
        colors="k",
        linewidths=0.5,
        levels=np.arange(0, 50, 10),
    )
    axs[row, 1].clabel(contour_tan, inline=True, fmt="%1.1f", fontsize=12)
    #[ax.invert_yaxis() for ax in (axs[row, 0], axs[row, 1])]


    # Set color bar labels
    axs[row, 0].set_title(f"{case_name} | $V_r$ | {prefix}")
    axs[row, 1].set_title(f"{case_name} | $V_t$| {prefix}")

    [ax.set_xlabel("Distance from the storm center (km)") for ax in (axs[row, 0], axs[row, 1])]
    [ax.set_ylabel("Pressure (hPa)") for ax in (axs[row, 0], axs[row, 1])]

    if row == 0:
     cb_rad = plt.colorbar(rad, ax=axs[:, 0].ravel(), orientation='horizontal')
     cb_tan = plt.colorbar(tan, ax=axs[:, 1].ravel(), orientation='horizontal')
     cb_rad.set_label(f"Radial Wind ($V_r$)")
     cb_tan.set_label(f"Tangential Wind ($V_t$)")
     return (cb_rad, cb_tan)
    else:
     None



pre = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2001/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-2*nc'))
post = sorted(glob.glob(
    '/nas/rgroup/stela/akumar/WRF_Harvey_v4/WRF_Simulations_FNL/LULC_2017/WRFV4.5.2/test/em_real/wrfout_d02_2017-08-2*nc'))

for index in range(len(pre)):
 prefix = pre[index].split('/')[-1][11:24]
 print(prefix)
 fig, axs = plt.subplots(2, 2, figsize=(14, 9.5), sharex=True, sharey=True)
 rad_tang_pre = xr.open_dataset(pre[index]).sel(level=slice(1000, 800))
 rad_tang_post = xr.open_dataset(post[index]).sel(level=slice(1000, 800))

 rad_c, tan_c = plot_rad_tang(rad_tang_pre, case_name = 'LULC 2001', prefix = prefix, row = 0)
 plot_rad_tang(rad_tang_post, case_name = 'LULC 2017', prefix = prefix, row = 1)
 plt.subplots_adjust(left=0.09, right=0.95, wspace=0.12, hspace=0.25, bottom=0.17, top=0.9)

 rad_c.ax.set_position([0.15, 0.07, 0.3, 0.02])
 tan_c.ax.set_position([0.58, 0.07, 0.3, 0.02])
 axs[0, 1].invert_yaxis()
 plt.savefig(f'../../figures_draft02/radial_tangential/{prefix}.jpeg', dpi=400, bbox_inches='tight')
 plt.close()
#plt.show()





