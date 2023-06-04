import matplotlib.pyplot as plt
import glob

import tropycal.tracks as tracks
import src.wrf_track

plt.rcParams.update({"font.size": 14, "font.weight": "bold"})


basin = tracks.TrackDataset(basin="north_atlantic", source="ibtracs", include_btk=False)
harvey = basin.get_storm(("harvey", 2017))

home_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS/'
wrf_runs = home_folder + "WRFV4_mp10_cu05/"
wrfoutfile1 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

homeic24_folder = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRF_mov1_GFS_IC24/'
wrf_runs = homeic24_folder + "WRFV4_mp10_cu06/"
wrfoutfile2 = sorted(glob.glob(wrf_runs + "wrfout_d03*"))

wrf_runs = '/nas/rstor/akumar/USA/PhD/Objective01/Hurricane_Harvey/WRF_Harvey_V2/WRF_Simulations/WRFV4_phy01_no_moving_domain/WRF/test/em_real/'
wrfoutfile5 = sorted(glob.glob(wrf_runs + "wrfout_d02*"))


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




out3 = src.wrf_track.get_track_details(wrfoutfile3[::2])
#out3a = src.wrf_track.get_track_details(wrfoutfile3a[::3])
#out3a = src.wrf_track.get_track_details(wrfoutfile3a[::3])
out3aa = src.wrf_track.get_track_details(wrfoutfile3aa)
out3ab = src.wrf_track.get_track_details(wrfoutfile3ab)
out3ac = src.wrf_track.get_track_details(wrfoutfile3ac)
#out7 = src.wrf_track.get_track_details(wrfoutfile7[::3])
#out10a = src.wrf_track.get_track_details(wrfoutfile10a)
#out10b = src.wrf_track.get_track_details(wrfoutfile10b)
#out10c = src.wrf_track.get_track_details(wrfoutfile10c)

#out = (out3, out10a, out10b, out10c)
#out = (out3, out3a, out3b)
out = (out3, out3aa, out3ab, out3ac)
#out = (out3,out8, out9 )
#out = (out3, out6, out6a) # out7, out8)
#lab = ('WRFV4_mp10_cu06', 'WRFV4_mp10_cu05 1 dom', 'WRFV4_mp10_cu05 2doms', '', '')
lab = ('WRFV4_mp10_cu06', '2dom cu01', '2dom cu05', '2dom cu06')
#lab = ('WRFV4_mp10_cu06', 'WRFV4_mp10_cu06 1dom', 'WRFV4_mp10_cu06 2 dom', '18km', 'post')
#lab = ('18km 1 dom', '18km 2 dom', '18km 3 dom')
col = ("b", "g", 'r', 'c', 'm', 'y')

src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)
 

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
lab = ('pre', 'post', 'new23', 'new19')
col = ('r', 'b', 'c', 'g')

src.wrf_track.plot_track_intensity(harvey, out, labels=lab, colors=col)
src.wrf_track.plot_track(harvey, out, labels=lab, colors=col)

plt.show()



"""



