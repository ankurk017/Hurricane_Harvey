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
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs
import matplotlib.patches as patches

plt.rcParams.update(
    {"font.size": 16, "font.weight": "bold", "savefig.dpi": 300})


def plot_domain(dbb, margin=4):
    fig = plt.figure(figsize=(9, 7))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines(resolution="10m")


    # Plot state coastlines
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none'
    )
    ax.add_feature(states_provinces, edgecolor='gray')
    


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
#    "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WPS/geo_em.d0*nc"))
    #"/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post/WPS/geo_em.d0*.nc"))
    "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WPS/geo_em.d0*.nc"))
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


url = "https://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2022-042723.txt"
basin = tracks.TrackDataset(basin='north_atlantic', atlantic_url=url)

cmap = plt.cm.jet
bounds = np.array([0, 32, 64, 83, 96, 113, 137, 150])
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

d1bb = get_bb(xr.open_dataset(geog_files[0]))
#d1bb[0] = -105
d2bb = get_bb(xr.open_dataset(geog_files[1]))
d3bb = get_bb(xr.open_dataset(geog_files[2]))
#d3bb = [-94, -89, 19.5, 24]
#d1bb = [-100, -85, 23, 35]


def plot_hurr_track(ax, year=2017, name='harvey', linestyle='dashed', color='k'):
    hurdat = basin.get_storm((name, year))
    ax.plot(hurdat['lon'], hurdat['lat'], color, linewidth=2, linestyle=linestyle, label=f'{name.capitalize()} ({year})')
    scatter = ax.scatter(hurdat['lon'], hurdat['lat'], 70,
                         c=hurdat['vmax'], cmap=cmap, norm=norm)
    return scatter


#ax = plot_domain((d1bb, d2bb, d3bb), margin=5)
ax = plot_domain((d1bb, d2bb, d3bb), margin=1)
# PLOT TRACK
shapefile_path = "/rhome/akumar/Downloads/Houston/COH_ADMINISTRATIVE_BOUNDARY_-_MIL.shp"
reader = shpreader.Reader(shapefile_path)
geometries = reader.geometries()

#coast.plot_coast(ax)
for geometry in geometries:
    ax.add_geometries([geometry], ccrs.PlateCarree(),
                         facecolor='none', edgecolor='blue')


scatter = plot_hurr_track(ax, year=2005, name='RITA', linestyle='dashed')
scatter = plot_hurr_track(ax, year=2008, name='IKE', linestyle='dotted')
scatter = plot_hurr_track(ax, year=2017, name='HARVEY', linestyle='solid')

plt.legend()

cbar = plt.colorbar(scatter, ax=ax, shrink=0.8)
cbar.set_ticks((bounds[:-1] + bounds[1:]) / 2)
cbar.set_ticklabels(["TD", "TS", "Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5"])

plt.tight_layout()
plt.savefig('../figures_draft/WRF_simulations_domain.jpeg')
plt.show()
