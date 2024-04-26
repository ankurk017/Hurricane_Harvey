import matplotlib as mpl
import tropycal.tracks as tracks
from netCDF4 import Dataset
from wrf import getvar, latlon_coords
import numpy as np
import matplotlib.pyplot as plt
import glob
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import matplotlib.patches as patches

plt.rcParams.update(
    {"font.size": 14, "font.weight": "bold", "savefig.dpi": 300})


def plot_domain(dbb, margin=4):
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines(resolution="10m")

    def draw_box(bb):
        return patches.Rectangle(
            (bb[0], bb[2]),
            bb[1] - bb[0],
            bb[3] - bb[2],
            linewidth=2,
            edgecolor="r",
            facecolor="none",
            alpha=1,
        )

    [ax.add_patch(draw_box(dbb)) for dbb in dbb]

    ax.xaxis.grid(True, which="major")
    ax.yaxis.grid(True, which="major")
    d1bb = dbb[0]
    ax.set_xlim(d1bb[0] - margin, d1bb[1] + margin)
    ax.set_ylim(d1bb[2] - margin, d1bb[3] + margin)

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1.5,
        color="k",
        alpha=0.45,
        linestyle="--",
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.right_labels = False
    return ax


geog_files = sorted(glob.glob(
    "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WPS/geo_em.d0*nc"))
A = xr.open_dataset(geog_files[0])


def get_bb(file_name):
    return np.array(
        (
            file_name["XLONG_M"].min().values,
            file_name["XLONG_M"].max().values,
            file_name["XLAT_M"].min().values,
            file_name["XLAT_M"].max().values,
        )
    )


# from wrfout
wrfoutfile = sorted(glob.glob(
    '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/wrfout_d01_2017-*'))


slp = getvar(Dataset(wrfoutfile[0]), "slp")
wrf_lat, wrf_lon = latlon_coords(slp)


basin = tracks.TrackDataset(basin='north_atlantic',
#                            source='hurdat', include_btk=False)
                            source='ibtracs', include_btk=False)

cmap = plt.cm.jet
bounds = [0, 32, 64, 83, 96, 113, 137, 150]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

d1bb = get_bb(xr.open_dataset(geog_files[0]))
d1bb[0] = -105
d2bb = get_bb(xr.open_dataset(geog_files[1]))
d3bb = [-94, -89, 19.5, 24]


def plot_hurr_track(ax, year=2017, name='harvey'):
    hurdat = basin.get_storm((name, year))
    ax.plot(hurdat['lon'], hurdat['lat'], "k")
    scatter = ax.scatter(hurdat['lon'], hurdat['lat'],
                         c=hurdat['vmax'], cmap=cmap, norm=norm)
    return scatter


ax = plot_domain((d1bb, d2bb, d3bb), margin=5)
# ax  =plot_domain((d1bb,), margin=5)

scatter = plot_hurr_track(ax, year=2017, name='harvey')
scatter = plot_hurr_track(ax, year=2021, name='ida')
scatter = plot_hurr_track(ax, year=2020, name='laura')
scatter = plot_hurr_track(ax, year=2020, name='delta')
# Landfall over Bahamas
scatter = plot_hurr_track(ax, year=2019, name='dorian')
scatter = plot_hurr_track(ax, year=2018, name='florence')
scatter = plot_hurr_track(ax, year=2018, name='michael')
scatter = plot_hurr_track(ax, year=2017, name='irma')
scatter = plot_hurr_track(ax, year=2016, name='matthew')

cbar = plt.colorbar(scatter, ax=ax)
cbar.ax.set_ylabel('10 m sustained Wind Speed (knots)')
plt.tight_layout()
#plt.savefig('../figures/Track_all_hurricanes.jpeg')
plt.show()
