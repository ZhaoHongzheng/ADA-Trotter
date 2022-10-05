# Plot the dynamics by ADA-Trotter and compare with exact quenched results
import numpy as np # generic math functions
import matplotlib
from colour import Color
from matplotlib import pyplot as plt
matplotlib.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1,2,figsize=(15,5))
for axs_j in axs:
    for axis in ['top', 'bottom', 'left', 'right']:
        axs_j.spines[axis].set_linewidth(3)

markersize = 12
linewidth = 4
interval=1
T_list = [1]
range_flag = 'short'
alpha = 3
GError=0.0
GError_random= 0.
J_zz_coupling = -1
h_x=-2
h_z = 0.2
theta = -np.pi/4
delta_x = 0
realization = 0
epsilon1 = 0.03
epsilon2_list = [0.1]
L_list = [24]
alpha_list = [1]

i=-1
j = 0
obs_str = 'H'
face_color_list = ['lightgreen','lightblue', 'white', 'k']
edge_color_list = ['darkgreen','mediumblue', 'k']
maker_list = ['-o', '-X','-d', '-', '-']
white = Color("red")
colors_blue = list(white.range_to(Color("blue"),3))
colors_red = list(white.range_to(Color("red"),len(epsilon2_list)+1))
f=-1
step_flag = '_N500'

for epsilon2 in epsilon2_list:
    j=j+1
    for L in L_list:
        file_name = 'data/AdaptiveBisection_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_L{L:.0f}_1epsilon{epsilon1:.4f}' \
                    '_2epsilon{epsilon2:.4f}_GError{GError:.3f}_RError{RError:.3f}_theta{theta:.2f}' \
                    '_deltax{delta_x:.2f}_r{r:.0f}'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z, L=L,
                                                           epsilon1=epsilon1
                                                           , epsilon2=epsilon2, GError=GError,
                                                           RError=GError_random, theta=theta,
                                                           delta_x=delta_x, r=realization) \
                    + range_flag + '_alpha' + str(alpha) + step_flag + '.npy'
        print(file_name)
        results_dic = np.load(file_name, allow_pickle=True)
        T_Flo = results_dic.item().get('T_Flo')
        obs_dic_RMD = results_dic.item().get('obs_dic_RMD')
        obs_dic_Static = results_dic.item().get('obs_dic_Static')
        obs_dic_eff = results_dic.item().get('obs_dic_eff')
        time_RMD_list = results_dic.item().get('time_RMD')
        time_Static_list = results_dic.item().get('time_Static')
        axs[0].plot(time_RMD_list, obs_dic_RMD['Mx'], '-o', markersize=markersize, linewidth=linewidth, zorder=20,
                    mfc='gold',
                    c='firebrick', alpha=1)
        axs[1].plot(time_RMD_list[0:-1],(results_dic.item().get('str_len_list')), '-', markersize=markersize, linewidth=3, zorder=20,
                    mfc='gold',
                    c='firebrick', alpha=1)

file_name_quench = 'data/Quench_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_T{T:.3f}_L{L:.0f}_GError{GError:.3f}_RError{RError:.3f}' \
                                 '_theta{theta:.2f}_'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z,
                                                             T=T_Flo, L=L, GError=GError, RError=GError_random,
                                                             theta=theta) + range_flag + '_alpha' + str(alpha) + step_flag +'_long.npy'
results_dic_Flo = np.load(file_name_quench, allow_pickle=True)
time_Flo_list = results_dic_Flo.item().get('time_Flo')
obs_dic_eff_Flo = results_dic_Flo.item().get('obs_dic_eff_Flo')
Floquet_step = T_Flo
linewidth = 4
axs[0].plot(time_Flo_list, obs_dic_eff_Flo['Mx'],'-', linewidth=linewidth, c='k',
                       alpha=1, zorder=5)

for epsilon2 in []:
    j=j+1
    for L in [24]:
        file_name = 'data/AdaptiveBisection_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_L{L:.0f}_1epsilon{epsilon1:.4f}' \
                    '_2epsilon{epsilon2:.4f}_GError{GError:.3f}_RError{RError:.3f}_theta{theta:.2f}' \
                    '_deltax{delta_x:.2f}_r{r:.0f}'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z, L=L,
                                                           epsilon1=epsilon1
                                                           , epsilon2=epsilon2, GError=GError,
                                                           RError=GError_random, theta=theta,
                                                           delta_x=delta_x, r=realization) \
                    + range_flag + '_alpha' + str(alpha) + '_N30' + '.npy'
        print(file_name)
        results_dic = np.load(file_name, allow_pickle=True)
        T_Flo = results_dic.item().get('T_Flo')
        obs_dic_RMD = results_dic.item().get('obs_dic_RMD')
        obs_dic_Static = results_dic.item().get('obs_dic_Static')
        obs_dic_eff = results_dic.item().get('obs_dic_eff')
        time_RMD_list = results_dic.item().get('time_RMD')
        time_Static_list = results_dic.item().get('time_Static')
        axs[0].plot(time_RMD_list, obs_dic_RMD['Mx'], '-D', markersize=markersize, linewidth=linewidth, zorder=2,
                    mfc='lightblue',
                    c='mediumblue', alpha=1)
        axs[1].plot(time_RMD_list[0:-1],(results_dic.item().get('str_len_list')), '-', markersize=markersize, linewidth=3, zorder=20,
                    mfc='gold',
                    c='firebrick', alpha=1)
axes = plt.gca()
plt.tick_params(labelsize=25,direction ='in', pad=7)
labelsize = 20
for axs_j in axs:
    axs_j.minorticks_on()
    axs_j.tick_params(which='major', length=15, width=2, direction='in', colors='k', labelsize=labelsize)
    axs_j.tick_params(which='minor', length=8, width=2, direction='in', colors='k', labelsize=labelsize)
    leg3 = axs_j.legend(prop={'size': 20}, frameon=False)
    leg3.get_frame().set_linewidth(1.5)
plt.show()
