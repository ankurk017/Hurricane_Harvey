import matplotlib.pyplot as plt
import glob
from src.wrf_src import wrf_assign_coords

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home = '/nas/rstor/xchen/class-project/Laura_WRF/'

wrf_pbl00 = sorted(glob.glob(home + f'/wrfout_d02*'))
#wrf_pbl00 = sorted(glob.glob(home + f'/WRF_FNL_post/WRF/test/em_real/wrfout_d02_2017-*'))
#wrf_pbl01 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/WRF_Building_Heights/WRF_Simulations/WRF/test/em_real/wrfout_d02_2017-*'))
#wrf_pbl02 = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/WRF_Building_Heights/WRF_Simulations/WRF_MYNN/test/em_real/wrfout_d02_2017-*'))

out_pbl1 = src.wrf_track.get_track_details(wrf_pbl00)
#out_pbl2 = src.wrf_track.get_track_details(wrf_pbl01)
#out_pbl3 = src.wrf_track.get_track_details(wrf_pbl02)


#out = (out_pbl1, out_pbl2, out_pbl3, )
out = (out_pbl1,  )

lab = ['2017',   ]
col = ("b", 'r','g', 'm')
titles = 'BH'


intensity_ax = src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig('../figures_paper/BH/pbl_track_pre_post.jpeg', dpi=400)
#intensity_ax.set_xlim(())

#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
#track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col, extent=[-103, -80, 16, 32])
track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col, extent=[-98.54, -91.72, 24.16, 32])
#plt.savefig('../figures_paper/BH/pbl_intensity_pre_post.jpeg', dpi=400)

#track_ax.set_xlim((-105, -85))
#track_ax.set_ylim((20, 45))
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
src.wrf_track.calculate_error(out, harvey)
pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures_paper/BH/{titles}_track_stats.csv')

plt.show()




