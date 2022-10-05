# Plot the long time error induced by ADA-Trotter, for different energy variance bounds
import numpy as np # generic math functions
import matplotlib
from colour import Color
from matplotlib import pyplot as plt
matplotlib.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1,2,figsize=(9,7))
for axs_j in axs:
    for axis in ['top', 'bottom', 'left', 'right']:
        axs_j.spines[axis].set_linewidth(2.5)

markersize = 12
linewidth=3
face_color_list = ['lightblue','lightgreen','gold', 'white', 'k']
edge_color_list = ['mediumblue','darkgreen','firebrick', 'k']

interval=1
T_list = [1]
range_flag = 'short'
alpha = 3
GError=0.0
GError_random= 0.
h_x=1
J_zz_coupling = 1
h_z = 0.3
delta_x = 0
realization = 0
d1_list = [0.001]
d2_list = [0,0.1,0.5,1,2,3,4,5,6,8,10]
theta_list = [np.pi/30,np.pi/8,np.pi/4,np.pi/5]
d1_index=0
theta_index = 0
theta = theta_list[theta_index]
L_list = [16,18,20]
alpha_list = [1,1,1]

i=-1
d2_index_list = [1,2,3,4,5,6,7,8]
white = Color("white")
colors_blue = list(white.range_to(Color("red"),len(L_list)+1))
j = -1

colors_list = colors_blue
for L in L_list:
    j += 1
    X_error_list = []
    Z_error_list = []
    epsilon2_list = []
    sd_z_list = []
    sd_x_list = []
    if L == 20:
        d2_index_list = [1, 2, 3, 4, 5,6]

    for d2_index in d2_index_list:
        epsilon1 = d1_list[d1_index]
        epsilon2 = d2_list[d2_index]
        try:
            file_name = 'data/Adaptive_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_L{L:.0f}' \
                        '_d1index{d1_index:.0f}_d2index{d2_index:.0f}_thetaindex{theta_index:.0f}' \
                        '.npy'.format(J=J_zz_coupling, h_x=h_x, h_z=h_z, L=L, d1_index=d1_index, d2_index=d2_index,
                                      theta_index=theta_index)

            results_dic = np.load(file_name, allow_pickle=True)
            obs_dic_RMD = results_dic.item().get('obs_dic_RMD')
            obs_dic_Static = results_dic.item().get('obs_dic_Static')
            obs_dic_eff = results_dic.item().get('obs_dic_eff')
            time_RMD_list = results_dic.item().get('time_RMD')
            time_Entropy_list = results_dic.item().get('time_Entropy')
            time_Static_list = results_dic.item().get('time_Static')
            epsilon2_list.append(epsilon2/L)

            sample_points = 500
            shift_z =  np.mean(obs_dic_eff['Mz'][-sample_points::])
            shift_x =  np.mean(obs_dic_eff['Mx'][-sample_points::])
            Z_error_list.append( np.abs(np.mean(obs_dic_RMD['Mz'][-sample_points::])-shift_z ))
            X_error_list.append( np.abs(np.mean(obs_dic_RMD['Mx'][-sample_points::])-shift_x ))
            sd_z_list.append(np.std(obs_dic_RMD['Mz'][-sample_points::]))
            sd_x_list.append(np.std(obs_dic_RMD['Mx'][-sample_points::]))

        except:
            continue
    axs[0].errorbar(x=epsilon2_list, y=Z_error_list,color = edge_color_list[j],mfc=face_color_list[j],
                       yerr=sd_z_list, capsize=0,linewidth=linewidth, marker='o',ms=markersize,alpha=alpha_list[j],
                       zorder=20)
    axs[1].errorbar(x=epsilon2_list, y=X_error_list,color = edge_color_list[j],mfc=face_color_list[j],
                    yerr=sd_x_list, capsize=0,linewidth=linewidth, marker='o',ms=markersize,alpha=alpha_list[j],
                    zorder=20)

axes = plt.gca()
plt.tick_params(labelsize=25,direction ='in', pad=7)

axs[0].set_xlabel("$d_2/L$", fontsize=25)
axs[1].set_xlabel("$d_2/L$", fontsize=25)
axs[0].set_ylabel("${O}_{\mathrm{AE}}-{O}_{\mathrm{DE}}$", fontsize=25)
labelsize=15
for axs_j in axs:
    axs_j.minorticks_on()
    axs_j.tick_params(which='major', length=15, width=2, direction='in', colors='k', labelsize=labelsize)
    axs_j.tick_params(which='minor', length=8, width=2, direction='in', colors='k', labelsize=labelsize)
    leg3 = axs_j.legend(prop={'size': 20}, frameon=False)
    leg3.get_frame().set_linewidth(1.5)

fig.suptitle('$J_z={Jz:.1f},h_x={hx:.1f}, h_z={h_z:0.1f}, d_1={epsilon1:.3f}$' \
 .format(Jz=J_zz_coupling,hx=h_x,  h_z=h_z,T_F=results_dic.item().get('T_Flo'), L=L,epsilon1=epsilon1), fontsize=25)

plt.show()
