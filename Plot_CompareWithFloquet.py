# Compare results obtained by using ADA-Trotter or Floquet dynamics.
import numpy as np # generic math functions
import matplotlib
from colour import Color
from matplotlib import pyplot as plt
matplotlib.rcParams['text.usetex'] = True
fig, axs = plt.subplots(2,2,figsize=(7,7))

for axs_i in axs:
    for axs_j in axs_i:
        for axis in ['top', 'bottom', 'left', 'right']:
            axs_j.spines[axis].set_linewidth(2.5)

def plot_bare_observable_simple_CompareFloquet(obs_all={},obs_eff={},str_len=0,time=0,time_Entropy=0
                                ,L=1,epsilon1=0,axs=0,edge_color='k',face_color='k',alpha=1,operator_str='Mx',marker='-o'):
    markersize=9
    Distance_list = []
    for t_i in range(len(time)):
        Distance_list.append((obs_all[operator_str][t_i] ))
    axs[0, 0].plot(time, Distance_list,marker,markersize=markersize, linewidth=3, zorder=20,mfc=face_color,
                   c=edge_color, alpha=alpha)
    axs[1, 0].plot(time,np.array(obs_all['H1']) / L,marker, markersize=markersize,mfc=face_color,
                   c=edge_color,
                   linewidth=3, alpha=alpha,zorder=10)

    Energy_variance = [(obs_all['H2'][i] - (obs_all['H1'][i]) ** 2) / L for i in range(len(time))]
    axs[1, 1].plot(time, Energy_variance, marker, markersize=markersize, linewidth=3,mfc=face_color,
                   c=edge_color,
    zorder = 10,alpha=alpha)

markersize = 9

interval=1
T_list = [1]
range_flag = 'short'
alpha = 3
str_option = ['+--+','-++-']
# str_option = ['+--+']
GError=0.0
GError_random= 0.0
h_x=-1.7
J_zz_coupling = -1
h_z = 0.5
theta = np.pi/8
delta_x = 0
realization = 0
epsilon1 = 0.03
epsilon2_list = [1]
L_list = [24]
alpha_list = [1]

i=-1
j = 0
obs_str = 'H'
face_color_list = ['lightblue', 'lightgreen','white', 'k']
edge_color_list = ['mediumblue','darkgreen', 'k']
maker_list = ['-o', '-X','-d', '-', '-']
white = Color("red")
colors_blue = list(white.range_to(Color("blue"),3))
colors_red = list(white.range_to(Color("red"),len(epsilon2_list)+1))
f=-1
step_flag = '_N30'
for epsilon2 in epsilon2_list:
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

        obs_dic_RMD = results_dic.item().get('obs_dic_RMD')
        obs_dic_Static = results_dic.item().get('obs_dic_Static')
        obs_dic_eff = results_dic.item().get('obs_dic_eff')
        time_RMD_list = results_dic.item().get('time_RMD')
        time_Static_list = results_dic.item().get('time_Static')

        for combo in [[time_RMD_list, obs_dic_RMD, obs_dic_eff],
                      ]:
            plot_bare_observable_simple_CompareFloquet(obs_all=combo[1], obs_eff=combo[2], time=combo[0]
                                                       , L=L, epsilon1=epsilon1, axs=axs, edge_color='firebrick',
                                                       face_color='gold', alpha=alpha_list[0])

            axs[0, 1].plot(time_RMD_list[1::], (results_dic.item().get('str_len_list')), '-o',
                           markersize=markersize, c='firebrick'
                           , linewidth=3, zorder=20, alpha=1, mfc='gold'
                           )
            axs[1, 0].plot(time_RMD_list, np.array(obs_dic_eff['H1']) / L, markersize=markersize, c='k',
                           linewidth=3, alpha=1, zorder=0)

            Energy_variance = [(obs_dic_eff['H2'][i] - (obs_dic_eff['H1'][i]) ** 2) / L for i in
                               range(len(time_RMD_list))]
            axs[1, 1].plot(time_RMD_list, Energy_variance, markersize=markersize, linewidth=3, color='k')

        T_flo = results_dic.item().get('T_Flo')
        for T_Flo in [0.08,T_flo]:
            if T_Flo == 0.08:
                marker = '-d'
            elif T_Flo == T_flo:
                marker = '->'
            i += 1
            file_name_Floshort = 'data/Floquet_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_T{T:.3f}_L{L:.0f}_GError{GError:.3f}_RError{RError:.3f}' \
                                 '_theta{theta:.2f}_'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z,
                                                             T=T_Flo, L=L, GError=GError, RError=GError_random,
                                                             theta=theta) + range_flag + '_alpha' + str(alpha) + step_flag +'.npy'

            print(file_name_Floshort)
            results_dic_Flo = np.load(file_name_Floshort, allow_pickle=True)
            obs_dic_Flo = results_dic_Flo.item().get('obs_dic_Flo')
            time_Flo_list = results_dic_Flo.item().get('time_Flo')
            obs_dic_eff_Flo = results_dic_Flo.item().get('obs_dic_eff_Flo')
            Floquet_step = T_Flo
            plot_bare_observable_simple_CompareFloquet(obs_all=obs_dic_Flo, obs_eff=obs_dic_eff_Flo, time= time_Flo_list
                                                       , L=L, epsilon1=epsilon1, axs=axs, edge_color=edge_color_list[i],face_color = face_color_list[i],
                                                       alpha=alpha_list[0],marker=marker)
            axs[0, 1].plot(time_Flo_list[1::], [(Floquet_step)] * len(time_Flo_list[1::]), marker,
                           markersize=markersize,mfc=face_color_list[i],
                   c=edge_color_list[i],  linewidth=3, zorder=30, alpha=alpha_list[0],
                           )

file_name_quench = 'data/Quench_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_T{T:.3f}_L{L:.0f}_GError{GError:.3f}_RError{RError:.3f}' \
                                 '_theta{theta:.2f}_'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z,
                                                             T=T_flo, L=L, GError=GError, RError=GError_random,
                                                             theta=theta) + range_flag + '_alpha' + str(alpha) + step_flag +'.npy'
print(file_name_quench)
results_dic_Flo = np.load(file_name_quench, allow_pickle=True)
time_Flo_list = results_dic_Flo.item().get('time_Flo')
obs_dic_eff_Flo = results_dic_Flo.item().get('obs_dic_eff_Flo')
Floquet_step = T_Flo
axs[0, 0].plot(time_Flo_list, obs_dic_eff_Flo['Mx'], linewidth=3, c='k',
                       alpha=0.8, zorder=30)


file_name = 'data/AdaptiveBisection_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_L{L:.0f}_1epsilon{epsilon1:.4f}' \
                    '_2epsilon{epsilon2:.4f}_GError{GError:.3f}_RError{RError:.3f}_theta{theta:.2f}' \
                    '_deltax{delta_x:.2f}_r{r:.0f}'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z, L=L,
                                                           epsilon1=epsilon1
                                                           , epsilon2=10, GError=GError,
                                                           RError=GError_random, theta=theta,
                                                           delta_x=delta_x, r=realization) \
                    + range_flag + '_alpha' + str(alpha) + step_flag + '.npy'
print(file_name)
results_dic = np.load(file_name, allow_pickle=True)
obs_dic_RMD = results_dic.item().get('obs_dic_RMD')
time_RMD_list = results_dic.item().get('time_RMD')

Energy_variance = [(obs_dic_RMD['H2'][i] - (obs_dic_RMD['H1'][i]) ** 2) / L for i in
                               range(len(time_RMD_list))]
axs[1,1].plot(time_RMD_list, Energy_variance, '--', markersize=markersize, linewidth=3, zorder=2,
                mfc='gold',
                c='firebrick', alpha=0.8)


axes = plt.gca()
plt.tick_params(labelsize=25,direction ='in', pad=7)

fontsize = 0
axs[1,0].set_xlabel("$t$", fontsize=fontsize)
axs[1,1].set_xlabel("$t$", fontsize=fontsize)
axs[0,0].set_ylabel("$\\Delta M_x$", fontsize=fontsize)
axs[1,1].set_ylabel("$(\\langle H^2\\rangle-\\langle H\\rangle^2)/{L}$", fontsize=fontsize)
axs[1,0].set_ylabel("$\\langle H \\rangle/L$", fontsize=fontsize)

fontsize = 25
for axs_i in axs:
    for axs_j in axs_i:
        axs_j.minorticks_on()
        axs_j.tick_params(which='major', length=15, width=2, direction='in', colors='k', labelsize=fontsize)
        axs_j.tick_params(which='minor', length=8, width=2, direction='in', colors='k', labelsize=fontsize)
        axs_j.set_xlim(0, 11)
        leg3 = axs_j.legend(prop={'size': 25}, frameon=False)
        leg3.get_frame().set_linewidth(1.5)

fig.suptitle('$J_z={Jz:.1f},h_x={hx:.1f}, h_z={h_z:0.1f}, L={L:.0f},d_1={epsilon1:.2f},d_2={epsilon2:.2f},$' \
 .format(Jz=J_zz_coupling,hx=h_x,  h_z=h_z,T_F=results_dic.item().get('T_Flo'), L=L,epsilon1=epsilon1,epsilon2=epsilon2), fontsize=fontsize)

plt.show()
