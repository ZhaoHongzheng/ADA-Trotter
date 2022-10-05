# Plot the violation of local Gauss's law for quantum link model.
import numpy as np # generic math functions
import matplotlib
from matplotlib import pyplot as plt
matplotlib.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1,2,figsize=(9,4))
for axs_j in axs:
    for axis in ['top', 'bottom', 'left', 'right']:
        axs_j.spines[axis].set_linewidth(2.5)

def plot_data(time_list,obs_dic, marker='-',linewidth=2,c='firebrick', mfc='gold',zorder=1,alpha=1,L=4,time_average_flag=False):
    Gauss1_error_list = []
    GaussV_error_list = []

    for t_i in range(len(time_list)):
        Gauss1_error = 0
        GaussV_error = 0
        for site in range(L):
            Gauss1_error = Gauss1_error + np.abs(
                obs_dic['Gauss1_' + str(site)][t_i] - obs_dic['Gauss1_' + str(site)][0]) / L
            GaussV_error = GaussV_error + np.abs(
                (obs_dic['Gauss2_' + str(site)][t_i] - obs_dic['Gauss1_' + str(site)][t_i] ** 2)
                - (obs_dic['Gauss2_' + str(site)][0] - obs_dic['Gauss1_' + str(site)][0] ** 2)) / L
        Gauss1_error_list.append(np.log10(Gauss1_error))
        GaussV_error_list.append(np.log10(GaussV_error))

    axs[0].plot(time_list, Gauss1_error_list, marker, linewidth=linewidth, c=c,mfc=mfc,alpha=alpha, zorder=zorder)
    axs[1].plot(time_list, GaussV_error_list, marker, linewidth=linewidth, c=c,mfc=mfc,alpha=alpha, zorder=zorder)


f=-1
str_max=300
A_Trotter = True
L_list = [6]
i=0
j = 0
obs_str = 'H'
face_color_list = ['gold','lightblue','lightgreen', 'white', 'k']
edge_color_list = ['firebrick','mediumblue','darkgreen', 'k']
maker_list = ['-o', '-X','-d', '-', '-']

J_0=0.5
if __name__ == "__main__":
    for gauge_breaking in [0.3]:
        for epsilon1, epsilon2 in [[0.1, 0.2]]:
            for epsilon_Gauss1, epsilon_GaussV in [[0.001, 0.003],[10,10]]:
                for L in L_list:
                    spin_S = 1
                    H_parm = {
                        'spin_S': spin_S,
                        'L': L,
                        'J': J_0 / ((spin_S * (spin_S + 1)) ** (1 / 2)),  # Matter-gauge couplings
                        'mu': 0.5,  # Gauge xx interaction
                        'k': 0.5,
                        'gauge_breaking_pert': gauge_breaking,  # Gauge breaking, z field on gauge sites
                        'gauge_preserving_pert': 0,  # Gauge x fields
                        'epsilon_Gauss1': epsilon_Gauss1,
                        'epsilon_GaussV': epsilon_GaussV,
                        'epsilon1': epsilon1,
                        'epsilon2': epsilon2,
                    }
                    if A_Trotter == True:
                        file_name = 'data/AdaptiveQLM_Sequential_J{J:.1f}_mu{mu:.2f}_k{k:.2f}_Gb{gauge_breaking_pert:.2f}_Gp{gauge_preserving_pert:.2f}_L{L:.0f}_1epsilon{epsilon1:.4f}' \
                                    '_2epsilon{epsilon2:.4f}_G1epsilon{epsilon_Gauss1:.4f}_GVepsilon{epsilon_GaussVariance:.4f}_Step{step:.0f}.npy'.format(
                            J=H_parm['J'],
                            mu=H_parm['mu'], L=L, k=H_parm['k'], epsilon_Gauss1=H_parm['epsilon_Gauss1'],
                            epsilon_GaussVariance=H_parm['epsilon_GaussV'],
                            epsilon1=H_parm['epsilon1'], epsilon2=H_parm['epsilon2'],
                            gauge_breaking_pert=H_parm['gauge_breaking_pert'],
                            gauge_preserving_pert=H_parm['gauge_preserving_pert'], step=str_max,
                        )

                        print(file_name)

                        results_dic = np.load(file_name, allow_pickle=True)
                        time_list = results_dic.item().get('time_RMD')
                        obs_dic = results_dic.item().get('obs_dic_RMD')
                        plot_data(time_list=time_list, obs_dic=obs_dic, marker='-', linewidth=3,c=edge_color_list[i], mfc=face_color_list[i],alpha=1, zorder=1,L=L)
                        i = i + 1
    c = 'k'
    mfc = 'k'
    alpha = 1
    zorder = 1

    axs[0].plot(time_list, [np.log10(0.001)]*len(time_list), '-', linewidth=2, c=c,
                   mfc=mfc, alpha=alpha, zorder=zorder)
    axs[1].plot(time_list, [np.log10(0.003)] * len(time_list), '-', linewidth=2, c=c,
                   mfc=mfc, alpha=alpha, zorder=zorder)


    axes = plt.gca()
    plt.tick_params(labelsize=25,direction ='in', pad=7)

    for axs_j in axs:
        axs_j.set_ylim(-3.8, -0.5)
        axs_j.minorticks_on()
        axs_j.tick_params(which='major', length=15, width=2, direction='in', colors='k', labelsize=25)
        axs_j.tick_params(which='minor', length=8, width=2, direction='in', colors='k', labelsize=25)
        axs_j.set_xlim(-1, 35)
        leg3 = axs_j.legend(prop={'size': 20}, frameon=False)
        leg3.get_frame().set_linewidth(1.5)

plt.show()
