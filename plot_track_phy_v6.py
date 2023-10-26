import matplotlib.pyplot as plt
import glob

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_FNL/'


home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_2dom/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))

wrfoutfile_ucm = sorted(glob.glob(home_2512 + f'/post_UCM/WRF//test/em_real/wrfout_d02_2017-*'))
wrfoutfile_urb = sorted(glob.glob(home_2512 + f'post_urban/WRF//test/em_real/wrfout_d02_2017-*'))

wrfoutfile_2512_def = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/def/WRF/test/em_real/wrfout_d02*'))
wrfoutfile_2512_pre = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/pre/WRF/test/em_real/wrfout_d02*'))
wrfoutfile_2512_post = sorted(glob.glob('/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2512/post/WRF_cntl/test/em_real/wrfout_d02*'))

out_wrfoutfile_pre_2612 = src.wrf_track.get_track_details(wrfoutfile_pre[::5])

#out_wrfoutfile_post = src.wrf_track.get_track_details(wrfoutfile_post[:50])
#out_wrfoutfile_ucm = src.wrf_track.get_track_details(wrfoutfile_ucm)
#out_wrfoutfile_urb = src.wrf_track.get_track_details(wrfoutfile_urb)
#out_wrfoutfile_def = src.wrf_track.get_track_details(wrfoutfile_2512_def[::3])
out_wrfoutfile_pre = src.wrf_track.get_track_details(wrfoutfile_2512_pre[::3])
out_wrfoutfile_pst = src.wrf_track.get_track_details(wrfoutfile_2512_post[::3])

#out = (out_wrfoutfile_pre, out_wrfoutfile_post)
#out = (out_wrfoutfile_post, out_wrfoutfile_ucm, out_wrfoutfile_urb)
#out = (out_wrfoutfile_post, out_wrfoutfile_urb)
out = (out_wrfoutfile_pre_2612, out_wrfoutfile_pre, out_wrfoutfile_pst)

lab = ('post', 'UCM', 'Urban', 'post 3km')
lab = ('2612', '2512_pre', '2512_post')

col = ("r", "g", 'b', 'm', 'c', 'y', 'olive', 'tomato', 'springgreen')
titles = 'set03'


src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
plt.savefig(f'../figures/{titles}_plot_intensity_2512.jpeg')
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
plt.savefig(f'../figures/{titles}_plot_track_2512.jpeg')
 
src.wrf_track.calculate_error(out, harvey)
pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures/{titles}_track_stats_2512.csv')

plt.show()

















