# Plot the long time error induced by ADA-Trotter, for different energy bounds
import numpy as np # generic math functions
import matplotlib
from colour import Color
from matplotlib import pyplot as plt
matplotlib.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1,2,figsize=(11,5))

for axs_j in axs:
    for axis in ['top', 'bottom', 'left', 'right']:
        axs_j.spines[axis].set_linewidth(2.5)

markersize = 12
linewidth=3

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
d1_list =[0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
d2_list = [0,0.1,0.5,1,2,3,4,5,6,8,10]
theta_list = [np.pi/30,np.pi/8,np.pi/4,np.pi/5]
d2_index=2
theta_index = 3
theta = theta_list[theta_index]
L_list = [16,18,20]
alpha_list = [1,1,0.7]

i=-1
d1_index_list = [1,2,3,4,5,6,7,8,9]
white = Color("white")
colors_blue = list(white.range_to(Color("red"),len(L_list)+1))

face_color_list = ['lightgreen','gold','lightblue', 'white', 'k']
edge_color_list = ['darkgreen','firebrick','mediumblue', 'k']

j = -1

colors_list = colors_blue
for L in L_list:
    j += 1
    X_error_list = []
    XX_error_list = []
    Z_error_list = []
    epsilon1_list = []
    epsilon2_list = []
    sd_z_list = []
    sd_x_list = []
    sd_xx_list = []

    for d1_index in d1_index_list:
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
            sample_points = 2000
            shift_z =  np.mean(obs_dic_eff['Mz'][-sample_points::])
            shift_xx = np.mean(obs_dic_eff['xx'][-sample_points::])
            shift_x =  np.mean(obs_dic_eff['Mx'][-sample_points::])

            Z_error_list.append( np.abs(np.mean(obs_dic_RMD['Mz'][-sample_points::])-shift_z ))
            X_error_list.append( np.abs(np.mean(obs_dic_RMD['Mx'][-sample_points::])-shift_x ))
            XX_error_list.append( np.abs(np.mean(obs_dic_RMD['xx'][-sample_points::])-shift_xx ))

            sd_z_list.append(np.std(obs_dic_RMD['Mz'][-sample_points::]))
            sd_x_list.append(np.std(obs_dic_RMD['Mx'][-sample_points::]))
            sd_xx_list.append(np.std(obs_dic_RMD['xx'][-sample_points::]))
            epsilon1_list.append(epsilon1)
        except:
            continue

    axs[0].errorbar(x=epsilon1_list, y=Z_error_list,color = edge_color_list[j],mfc=face_color_list[j],
                       yerr=sd_z_list, capsize=0,linewidth=linewidth, marker='o',ms=markersize,
                       zorder=20)
    axs[1].errorbar(x=epsilon1_list, y=X_error_list,color = edge_color_list[j],mfc=face_color_list[j],
                    yerr=sd_x_list, capsize=0,linewidth=linewidth, marker='o',ms=markersize,
                    zorder=20)

axes = plt.gca()
plt.tick_params(labelsize=25,direction ='in', pad=7)

axs[0].set_xlabel("$d_1$", fontsize=25)
axs[1].set_xlabel("$d_1$", fontsize=25)
axs[0].set_ylabel("${O}_{\mathrm{AE}}-{O}_{\mathrm{DE}}$", fontsize=25)


L=20
X_micro_error_list = []
XX_micro_error_list = []
Z_micro_error_list = []
epsilon1_list = []
file_name = 'data/Spectrum_J{J:.1f}_hx{h_x:.2f}_hz{h_z:0.2f}_L{L:.0f}_theta'.format(J=J_zz_coupling, h_x=h_x,
     h_z=h_z, L=L) + str(theta) + '.npy'
results_dic = np.load(file_name, allow_pickle=True)
micro_prediction_z = results_dic.item().get('Mz_micro')
micro_prediction_x = results_dic.item().get('Mx_micro')
micro_prediction_xx = results_dic.item().get('xx_micro')
for d1_index in d1_index_list:
    epsilon1 = d1_list[d1_index]
    epsilon2 = d2_list[d2_index]
    epsilon1_list.append(epsilon1)
    X_micro_error_list.append(np.abs(micro_prediction_x[str(epsilon1)] - micro_prediction_x[str(0)]))
    Z_micro_error_list.append(np.abs(micro_prediction_z[str(epsilon1)] - micro_prediction_z[str(0)]))
    XX_micro_error_list.append(np.abs(micro_prediction_xx[str(epsilon1)] - micro_prediction_xx[str(0)]))

axs[0].plot(epsilon1_list, Z_micro_error_list, 'D-', color='grey',mec='k',
                               markersize=markersize, linewidth=linewidth, zorder=0, alpha=1
                               )
axs[1].plot(epsilon1_list, X_micro_error_list, 'D-', color='grey', mec='k',
            markersize=markersize, linewidth=linewidth, zorder=0, alpha=1
            )
labelsize = 15
for axs_j in axs:
    axs_j.minorticks_on()
    axs_j.tick_params(which='major', length=15, width=2, direction='in', colors='k', labelsize=labelsize)
    axs_j.tick_params(which='minor', length=8, width=2, direction='in', colors='k', labelsize=labelsize)
    leg3 = axs_j.legend(prop={'size': 20}, frameon=False)
    leg3.get_frame().set_linewidth(1.5)

fig.suptitle('$J_z={Jz:.1f},h_x={hx:.1f}, h_z={h_z:0.1f}, d_2={epsilon2:.3f}$' \
 .format(Jz=J_zz_coupling,hx=h_x,  h_z=h_z,T_F=results_dic.item().get('T_Flo'), L=L,epsilon2=epsilon2), fontsize=25)

plt.show()
