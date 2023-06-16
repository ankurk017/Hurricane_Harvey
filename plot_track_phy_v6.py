import matplotlib.pyplot as plt
import glob

import tropycal.tracks as tracks
import src.wrf_track
import pandas as pd
plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_FNL/'

wrfoutfile1a = sorted(glob.glob(home_folder + "WRFV4_mp06_cu01/" + "wrfout_d03*"))[::1]
wrfoutfile1b = sorted(glob.glob(home_folder + "WRFV4_mp06_cu05/" + "wrfout_d03*"))[::1]
wrfoutfile1c = sorted(glob.glob(home_folder + "WRFV4_mp06_cu06/" + "wrfout_d03*"))[::1]
wrfoutfile1d = sorted(glob.glob(home_folder + "WRFV4_mp08_cu01/" + "wrfout_d03*"))[::1]
wrfoutfile1e = sorted(glob.glob(home_folder + "WRFV4_mp08_cu05/" + "wrfout_d03*"))[::1]
wrfoutfile1f = sorted(glob.glob(home_folder + "WRFV4_mp08_cu06/" + "wrfout_d03*"))[::1]
wrfoutfile1g = sorted(glob.glob(home_folder + "WRFV4_mp10_cu01/" + "wrfout_d03*"))[::1]
wrfoutfile1h = sorted(glob.glob(home_folder + "WRFV4_mp10_cu05/" + "wrfout_d03*"))[::1]
wrfoutfile1i = sorted(glob.glob(home_folder + "WRFV4_mp10_cu06/" + "wrfout_d03*"))[::1]

homeic24_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC24/'
wrf_runs = homeic24_folder + "WRFV4_mp10_cu06/"
wrfoutfile2 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/'

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain/WRF/test/em_real/'
wrfoutfile5 = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06/'
wrfoutfile3 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06_1doms_BEM/'
wrfoutfile3a = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC/WRFV4_mp10_cu06_2doms_BEM/'
wrfoutfile3b = sorted(glob.glob(wrf_runs + "wrfout_d02*"))


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_v2/def/WRF_mp10_cu01/'
wrfoutfile3aa = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_v2/def/WRF_mp10_cu05/'
wrfoutfile3ab = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_v2/def/WRF_mp10_cu05_no_ocean_physics/'
wrfoutfile3ab_noocean = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_v2/def/WRF_mp10_cu06/'
wrfoutfile3ac = sorted(glob.glob(wrf_runs + "wrfout_d02*"))


wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain/WRF_2dom_BEP/test/em_real/'
wrfoutfile6 = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain/WRF/test/em_real/'
wrfoutfile6a = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_18km/WRFV4_18km_dom01/'
wrfoutfile7 = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain_pre-post/pre/WRF/test/em_real/'
wrfoutfile8 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain_pre-post/post/WRF/test/em_real/'
wrfoutfile9 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_18km/def/WRF_1dom/'
wrfoutfile10a = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_18km/def/WRF_2dom/'
wrfoutfile10b = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC25_12UTC_18km/def/WRF_3dom/'
wrfoutfile10c = sorted(glob.glob(wrf_runs + "wrfout_d03*"))


wrf_runs = home_2512 + '/WRF_mov1_GFS_IC25_12UTC_v2/def/WRF_mp10_cu05_no_ocean_physics/'
wrfoutfile11def = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrf_runs = home_2512 + '/WRF_mov1_GFS_IC25_12UTC_v2/pre/WRF_mp10_cu05_no_ocean_physics/'
wrfoutfile11pre = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrf_runs = home_2512 + '/WRF_mov1_GFS_IC25_12UTC_v2/post/WRF_mp10_cu05_no_ocean_physics/'
wrfoutfile11post = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrf_runs = home_2512 + '/WRF_mov1_GFS_IC26_00UTC/WRF_mp10_cu05_no_ocean_physics/'
wrfoutfile6_ic2600 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = home_2512 + 'WRFV45/WRF/test/em_real/'
wrfoutfile_wrf45 = sorted(glob.glob(wrf_runs + "wrfout_d01*"))
wrf_runs = home_2512 + 'WRFV45_FNL/WRF/test/em_real/old/'
wrf_runs = home_2512 + 'WRFV45_FNL/WRF/test/em_real/'
wrf_runs = home_2512 + 'WRFV45_FNL_2518/WRF/test/em_real/'
wrf_runs = home_2512 + 'WRFV45_FNL_2518/WRF_old/test/em_real/'
wrfoutfile_wrf45_fnl = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

wrf_runs = home_2512 + 'WRF_mov1_FNL_2512/WRFV4_mp10_cu06/'
wrfoutfile_wrf43_fnl = sorted(glob.glob(wrf_runs + "wrfout_d01*"))[::5]


wrf_runs = home_2512 + 'WRFV45_FNL_2518/WRF_old_2doms/test/em_real/'
wrfoutfile_wrf43_fnl_dom2 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))
wrfoutfile_wrf43_fnl_dom1 = sorted(glob.glob(wrf_runs + "wrfout_d01*"))

#wrf_runs = home_2512 + 'WRF_FNL_2612/WRF/test/em_real/'
#wrfoutfile_ic2612a = sorted(glob.glob(wrf_runs + "wrfout_d01*"))
#wrfoutfile_ic2612b = sorted(glob.glob(wrf_runs + "wrfout_d02*"))


home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V3/WRF_Simulations/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/WRF_HRRR_sim01_pre/WRF/test/em_real/wrfout_d01_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/WRF_HRRR_sim01_post/WRF/test/em_real/wrfout_d01_2017-*'))

home_2512 = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_FNL_2612/'
wrfoutfile_pre = sorted(glob.glob(home_2512 + f'/pre/WRF_2dom/test/em_real/wrfout_d02_2017-*'))
wrfoutfile_post = sorted(glob.glob(home_2512 + f'/post/WRF_2dom/test/em_real/wrfout_d02_2017-*'))


#out1a = src.wrf_track.get_track_details(wrfoutfile1a)
#out1b = src.wrf_track.get_track_details(wrfoutfile1b)
#out1c = src.wrf_track.get_track_details(wrfoutfile1c)
#out1d = src.wrf_track.get_track_details(wrfoutfile1d)
#out1e = src.wrf_track.get_track_details(wrfoutfile1e)
#out1f = src.wrf_track.get_track_details(wrfoutfile1f)
#out1g = src.wrf_track.get_track_details(wrfoutfile1g)
#out1h = src.wrf_track.get_track_details(wrfoutfile1h)
#out1i = src.wrf_track.get_track_details(wrfoutfile1i)

#out3 = src.wrf_track.get_track_details(wrfoutfile3[::2])
#out3a = src.wrf_track.get_track_details(wrfoutfile3a[::3])
#out3b = src.wrf_track.get_track_details(wrfoutfile3b[::3])
#out3aa = src.wrf_track.get_track_details(wrfoutfile3aa[::3])
#out3ab = src.wrf_track.get_track_details(wrfoutfile3ab[::3])
#out3ab_noocean = src.wrf_track.get_track_details(wrfoutfile3ab_noocean)
#out3ac = src.wrf_track.get_track_details(wrfoutfile3ac[::3])
#out5 = src.wrf_track.get_track_details(wrfoutfile5)
#out6 = src.wrf_track.get_track_details(wrfoutfile6[::3])
#out6a = src.wrf_track.get_track_details(wrfoutfile6a[::2])
#out7 = src.wrf_track.get_track_details(wrfoutfile7[::3])
#out10a = src.wrf_track.get_track_details(wrfoutfile10a)
#out10b = src.wrf_track.get_track_details(wrfoutfile10b)
#out10c = src.wrf_track.get_track_details(wrfoutfile10c)
#out11_def = src.wrf_track.get_track_details(wrfoutfile11def[::1])
#out11_pre = src.wrf_track.get_track_details(wrfoutfile11pre[::1])
#out11_post = src.wrf_track.get_track_details(wrfoutfile11post[::1])
#out12_ic26 = src.wrf_track.get_track_details(wrfoutfile6_ic2600)

#out2518pre = src.wrf_track.get_track_details(wrfoutfile_ic2518pre)
#out2518post = src.wrf_track.get_track_details(wrfoutfile_ic2518post)
#out_wrfoutfile_wrf45 = src.wrf_track.get_track_details(wrfoutfile_wrf45)
#out_wrfoutfile_wrf45_fnl = src.wrf_track.get_track_details(wrfoutfile_wrf45_fnl)
#out_wrfoutfile_wrf43_fnl = src.wrf_track.get_track_details(wrfoutfile_wrf43_fnl)
#out_wrfoutfile_wrf45_fnl_dom2 = src.wrf_track.get_track_details(wrfoutfile_wrf43_fnl_dom2)
#out_wrfoutfile_wrf45_fnl_dom1 = src.wrf_track.get_track_details(wrfoutfile_wrf43_fnl_dom1)
#out_wrfoutfile_ic2612a = src.wrf_track.get_track_details(wrfoutfile_ic2612a)
#out_wrfoutfile_ic2612b = src.wrf_track.get_track_details(wrfoutfile_ic2612b)

out_wrfoutfile_pre = src.wrf_track.get_track_details(wrfoutfile_pre)
out_wrfoutfile_post = src.wrf_track.get_track_details(wrfoutfile_post)
out = (out_wrfoutfile_pre, out_wrfoutfile_post)

#out = (out2518pre, out2518post)
#out = (out_wrfoutfile_wrf45, out_wrfoutfile_wrf45_fnl, out_wrfoutfile_wrf43_fnl)
#out = (out_wrfoutfile_wrf45_fnl, out_wrfoutfile_wrf43_fnl)
#out = (out_wrfoutfile_wrf45_fnl_dom1, out_wrfoutfile_wrf45_fnl_dom2)
#out = (out_wrfoutfile_ic2612a,  out_wrfoutfile_ic2612b)
#out = (out_wrfoutfile_ic2612a,  )
#out = (out7, out10a, out10b, out10c)
#out = (out3, out3a, out3b) # works like a charm  ########## YEAHHHHHHHH
#out = (out3, out3aa, out3ab, out3ac) # this is also good 
#out = (out3, out3aa, out3ab, out3ab_noocean) # this is also good 
#out = (out3, out11_def, out11_pre, out11_post, out12_ic26) # latest siulations - best pre and post
#out = (out11_def, out11_pre, out11_post) # latest siulations - best pre and post
#out = (out3,out8, out9 )
#out = (out3, out6, out6a) # out7, out8)
#out = (out5, out6, out6a) # out7, out8)
#out = (out1a, out1b, out1c, out1d, out1e, out1f, out1g, out1h, out1i,)
#out = (out5, ) # out7, out8)
#lab = ('WRFV4_mp10_cu06', 'WRFV4_mp10_cu05 1 dom', 'WRFV4_mp10_cu05 2doms', '', '')
#lab = ('WRFV4_mp10_cu06', '2dom cu01', '2dom cu05', '2dom cu05_nooncea')
#lab = ('3doms', '2 doms def', '2 doms pre', '2 doms post', 'IC2600')
lab = ('def 3km', 'pre 3km', 'post 3km')
lab = ('IC2518_GFS', 'IC2518_FNL', 'post 3km')
#lab = ('WRFV4_mp06_cu01', 'WRFV4_mp06_cu05', 'WRFV4_mp06_cu06', 'WRFV4_mp08_cu01', 'WRFV4_mp08_cu05', 'WRFV4_mp08_cu06', 'WRFV4_mp10_cu01', 'WRFV4_mp10_cu05', 'WRFV4_mp10_cu06',)
#lab = ('WRFV4_mp10_cu06', 'WRFV4_mp10_cu06 1dom', 'WRFV4_mp10_cu06 2 dom', '18km', 'post')
#lab = ('18km 1 dom', '18km 2 dom', '18km 3 dom')
col = ("r", "g", 'b', 'm', 'c', 'y', 'olive', 'tomato', 'springgreen')
titles = 'set03'



src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_intensity.jpeg')
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
#plt.savefig(f'../figures/{titles}_plot_track.jpeg')
 
#src.wrf_track.calculate_error(out, harvey)
#pd.DataFrame(src.wrf_track.calculate_error(out, harvey), index=lab[:len(out)]).to_csv(f'../figures/{titles}_track_stats.csv')

plt.show()

"""

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4/'
wrfoutfile1 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/post/WRF_9-3-51/WRFV4/'
wrfoutfile2 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4_test_may23/'
wrfoutfile3 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V1/pre/WRF_9-3-51/WRFV4_test_may19/'
wrfoutfile4 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))

out1 = src.wrf_track.get_track_details(wrfoutfile1[::3])
out2 = src.wrf_track.get_track_details(wrfoutfile2[::3])
out3 = src.wrf_track.get_track_details(wrfoutfile3[::3])
out4 = src.wrf_track.get_track_details(wrfoutfile4[::3])

out = (out1, out2, out3, out4)
lab = ('IC24_pre', 'IC24_post', 'IC24 cu05', 'IC24 cu06')
col = ('r', 'b', 'c', 'g')

src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)

plt.savefig('../figures/set01_plot_track.jpeg')
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
plt.savefig('../figures/set01_plot_intensity.jpeg')

plt.show()

"""
















