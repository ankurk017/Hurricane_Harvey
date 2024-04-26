import matplotlib.pyplot as plt
import glob
from src.wrf_src import wrf_assign_coords

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 16, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/"
wrf_pbl01 = sorted(glob.glob(f'{home}/wrfout_d02_2017-*'))[::6]

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF/test/em_real/"
wrf_pbl02 = sorted(glob.glob(f'{home}/wrfout_d02_2017-*'))[::6]

home = "/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF_vortex/test/em_real/"
wrf_pbl02 = sorted(glob.glob(f'{home}/wrfout_d02_2017-*'))[::6]


out_pbl1 = src.wrf_track.get_track_details(wrf_pbl01)
out_pbl2 = src.wrf_track.get_track_details(wrf_pbl02)


out = (out_pbl1, out_pbl2)

lab = ('LULC 2001', 'Vortex', )
col = ("b", 'r',)
titles = 'LULC'


intensity_ax = src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig('../figures_paper/track_pre_post.jpeg', dpi=300)
#intensity_ax.set_xlim(())

#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
track_ax = src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
#plt.savefig('../figures_paper/intensity_pre_post.jpeg', dpi=300)


#track_ax.set_xlim((-105, -85))
#track_ax.set_ylim((20, 45))
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
#src.wrf_track.calculate_error(out, harvey)
#pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures_paper/{titles}_track_stats.csv')

plt.show()



