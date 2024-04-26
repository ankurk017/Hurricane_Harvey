import matplotlib.pyplot as plt
import glob

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post_ensemble/'

wrf_plu1 = sorted(glob.glob(home_2512 + f'WRF_plu1/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_plu2 = sorted(glob.glob(home_2512 + f'WRF_plu2/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_plu3 = sorted(glob.glob(home_2512 + f'WRF_plu3/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_min1 = sorted(glob.glob(home_2512 + f'WRF_min1/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_min2 = sorted(glob.glob(home_2512 + f'WRF_min2/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_min3 = sorted(glob.glob(home_2512 + f'WRF_min3/test/em_real/wrfout_d02_2017-*'))[:80][::3]
wrf_cntl = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))[:80][::3]




out_cntl = src.wrf_track.get_track_details(wrf_cntl)
out_plu1 = src.wrf_track.get_track_details(wrf_plu1)
out_plu2 = src.wrf_track.get_track_details(wrf_plu2)
out_plu3 = src.wrf_track.get_track_details(wrf_plu3)
out_min1 = src.wrf_track.get_track_details(wrf_min1)
out_min2 = src.wrf_track.get_track_details(wrf_min2)
out_min3 = src.wrf_track.get_track_details(wrf_min3)

out = (out_cntl, out_plu1, out_plu2, out_plu3, out_min1, out_min2, out_min3)

lab = ('cntl', '+1%', '+2%', '+3%', '-1%', '-2%', '-3%')

col = ("r", "g", 'b', 'm', 'c', 'y', 'olive', 'tomato', 'springgreen')
titles = 'set03'


src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
#src.wrf_track.calculate_error(out, harvey)
#pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures/{titles}_track_stats.csv')

plt.show()

















