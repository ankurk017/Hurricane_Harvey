import matplotlib.pyplot as plt
import glob

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post_ensemble/'
wrf_cntl = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))[:80][::3]




out_cntl = src.wrf_track.get_track_details(wrf_cntl)


out = (out_cntl, )

lab = ('cntl', )

col = ("r", )
titles = 'set03'


src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
#src.wrf_track.calculate_error(out, harvey)
#pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures/{titles}_track_stats.csv')

plt.show()

















