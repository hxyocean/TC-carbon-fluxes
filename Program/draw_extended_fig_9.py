import numpy as np
import warnings
import time as Time
import netCDF4 as nc
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scienceplots
import scipy.io as scio


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    time_start = Time.time()
    matplotlib.style.use(['science', 'nature', 'no-latex'])
    matplotlib.rcParams['font.family'] = 'Verdana'


    lw = 2
    lsize = 17
    m_a = 0.8
    font_l = {'family': 'Verdana', 'weight': 'bold', 'size': lsize + 2, }
    c1 = ['#1762a1', '#01847F', 'gray', '#703274', '#c65127', '#e92921']


    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_17/pco2.mat')
    dpco2_all = np.squeeze(tem_data['dpco2_all'])
    var_F_SST = np.squeeze(tem_data['var_F_SST'])
    var_F_DIC = np.squeeze(tem_data['var_F_DIC'])
    var_F_TA = np.squeeze(tem_data['var_F_TA'])
    var_SST = np.squeeze(tem_data['var_SST'])
    var_DIC = np.squeeze(tem_data['var_DIC'])
    var_TA = np.squeeze(tem_data['var_TA'])
    v_pre_SST = np.squeeze(tem_data['v_pre_SST'])
    v_pre_DIC = np.squeeze(tem_data['v_pre_DIC'])
    v_pre_TA = np.squeeze(tem_data['v_pre_TA'])
    F_SST = np.squeeze(tem_data['F_SST'])
    F_DIC = np.squeeze(tem_data['F_DIC'])
    F_TA = np.squeeze(tem_data['F_TA'])

    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_17/F_line.mat')
    F_line = np.squeeze(tem_data['F_line'])
    T_line_kpp = np.squeeze(tem_data['T_line_kpp'])
    T_line_prt = np.squeeze(tem_data['T_line_prt'])
    T_line_pwp = np.squeeze(tem_data['T_line_pwp'])


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    fig = plt.figure(figsize=(20, 15), dpi=250)
    G = gridspec.GridSpec(3, 6, figure=fig)
    ax1 = plt.subplot(G[0, :])
    plt.title('a', x=-0.064, y=1.1, font=font_l, horizontalalignment='left', verticalalignment='top')
    plt.axis('off')
    ax = ax1.inset_axes([0., 0., 1, 1.])
    l1, = ax.plot(-F_line, T_line_prt, linewidth=lw, linestyle='-', alpha=m_a, color=c1[3])
    ax.plot(-F_line, T_line_prt, linewidth=lw, linestyle='', marker='o', markersize=8, alpha=0.5, color=c1[3])
    l2, = ax.plot(-F_line, T_line_kpp, linewidth=lw, linestyle='-', alpha=m_a, color=c1[1])
    ax.plot(-F_line, T_line_kpp, linewidth=lw, linestyle='', marker='o', markersize=8, alpha=0.5, color=c1[1])
    l3, = ax.plot(-F_line, T_line_pwp, linewidth=lw, linestyle='-', alpha=m_a, color=c1[4])
    ax.plot(-F_line, T_line_pwp, linewidth=lw, linestyle='', marker='o', markersize=8, alpha=0.5, color=c1[4])
    ax.tick_params(labelsize=lsize)
    ax.grid(linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
    ax.set_xlabel(r'$\Gamma_{T}\ (^{\circ} C\ m^{-1})$', fontsize=lsize - 1)
    ax.set_ylabel(r'$SST\ anomaly\ (^{\circ} C)$', fontsize=lsize - 1)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.002))
    ax.legend([l2, l3, l1], ['KPP', 'PWP', 'PRT'], fontsize=lsize, frameon=True, framealpha=0.6)
    ax.text(x=-0.0498, y=-1.44, s=r'$Strengthening$', fontsize=lsize, horizontalalignment='center', verticalalignment='top')
    ax.quiver(-0.0478, -1.49, -1, 0, color='k', angles='xy', scale=5, width=0.008, headwidth=2.5, headlength=5, alpha=0.8)
    ax.set_ylim([-1.52, -0.98])

    ax = plt.subplot(G[1, 0:3])
    l1, = plt.plot(v_pre_DIC + var_DIC, np.abs((dpco2_all[:, 0] - dpco2_all[:, 1]) / dpco2_all[:, 1]) * 100,
                   linestyle='-', linewidth=lw, alpha=m_a, color=c1[3])
    plt.tick_params(labelsize=lsize)
    plt.grid(linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
    plt.ylim([0, 2])
    plt.ylabel(r'$Error\ percentage\ (\%)$', fontsize=lsize)
    plt.text(1.03, -0.05, r'$DIC_{0}\ (\mu mol\ kg^{-1})$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[3], ha='left')
    ax.xaxis.set_major_locator(plt.MultipleLocator(100))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))

    ax1 = ax.twiny()
    plt.plot(v_pre_TA + var_TA, np.abs((dpco2_all[:, 2] - dpco2_all[:, 3]) / dpco2_all[:, 3]) * 100,
             linestyle='-', linewidth=lw, alpha=m_a, color=c1[1])
    plt.tick_params(labelsize=lsize)
    ax1.spines["bottom"].set_position(("outward", 25))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_label_position('bottom')
    plt.text(1.03, -0.15, r'$TA_{0}\ (\mu mol\ kg^{-1})$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[1], ha='left')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(100))

    ax2 = ax.twiny()
    plt.plot(v_pre_SST + var_SST, np.abs((dpco2_all[:, 4] - dpco2_all[:, 5]) / dpco2_all[:, 5]) * 100,
             linestyle='-', linewidth=lw, alpha=m_a, color=c1[4])
    plt.tick_params(labelsize=lsize)
    plt.text(1.03, 1.03, r'$T_{0}\ (^{\circ} C)$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[4], ha='left')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(2.5))

    ax.spines['top'].set_color(c1[4])
    ax.spines['bottom'].set_color(c1[3])
    ax1.spines['top'].set_color(c1[4])
    ax1.spines['bottom'].set_color(c1[1])
    ax2.spines['top'].set_color(c1[4])
    ax2.spines['bottom'].set_color(c1[3])
    ax.tick_params(axis='x', colors=c1[3])
    ax1.tick_params(axis='x', colors=c1[1])
    ax2.tick_params(axis='x', colors=c1[4])
    plt.title('b', x=-0.15, y=1.1, font=font_l, horizontalalignment='left',
              verticalalignment='top')


    ax = plt.subplot(G[1, 3:])
    l1, = plt.plot(-(F_DIC + var_F_DIC), np.abs((dpco2_all[:, 8] - dpco2_all[:, 9]) / dpco2_all[:, 9]) * 100,
                   linestyle='-', linewidth=lw, alpha=m_a, color=c1[3])
    plt.tick_params(labelsize=lsize)
    plt.grid(linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
    plt.ylim([0, 2])
    plt.ylabel(r'$Error\ percentage\ (\%)$', fontsize=lsize)
    plt.text(1.03, -0.05, r'$\Gamma_{DIC}\ (\mu mol\ kg^{-1}\ m^{-1})$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[3], ha='left')
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))

    ax1 = ax.twiny()
    plt.plot(-(F_TA + var_F_TA), np.abs((dpco2_all[:, 10] - dpco2_all[:, 11]) / dpco2_all[:, 11]) * 100,
             linestyle='-', linewidth=lw, alpha=m_a, color=c1[1])
    plt.tick_params(labelsize=lsize)
    ax1.spines["bottom"].set_position(("outward", 25))
    ax1.xaxis.set_ticks_position('bottom')
    ax1.xaxis.set_label_position('bottom')
    plt.text(1.03, -0.15, r'$\Gamma_{TA}\ (\mu mol\ kg^{-1}\ m^{-1})$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[1], ha='left')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(0.05))

    ax2 = ax.twiny()
    plt.plot(-(F_SST + var_F_SST), np.abs((dpco2_all[:, 12] - dpco2_all[:, 13]) / dpco2_all[:, 13]) * 100,
             linestyle='-', linewidth=lw, alpha=m_a, color=c1[4])
    plt.tick_params(labelsize=lsize)
    plt.text(1.06, 1.03, r'$\Gamma_{T}\ (^{\circ} C\ m^{-1})$', transform=ax.transAxes, fontsize=lsize - 2, color=c1[4], ha='left')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(0.005))

    ax.spines['top'].set_color(c1[4])
    ax.spines['bottom'].set_color(c1[3])
    ax1.spines['top'].set_color(c1[4])
    ax1.spines['bottom'].set_color(c1[1])
    ax2.spines['top'].set_color(c1[4])
    ax2.spines['bottom'].set_color(c1[3])
    ax.tick_params(axis='x', colors=c1[3])
    ax1.tick_params(axis='x', colors=c1[1])
    ax2.tick_params(axis='x', colors=c1[4])

    plt.title('c', x=-0.15, y=1.1, font=font_l, horizontalalignment='left',
              verticalalignment='top')


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    tem_x_label = [r'$DIC_{0}$', r'$TA_{0}$', r'$T_{0}$',
                   r'$\Gamma_{DIC}$', r'$\Gamma_{TA}$', r'$\Gamma_{T}$',]
    m_error_abs = []
    for i in range(8):
        tem_data = dpco2_all[:, i*2 + 1] - dpco2_all[:, i*2]
        m_error_abs.append(tem_data / dpco2_all[:, i*2 + 1] * 1e2)
    tem_x = np.linspace(1, 6, 6)
    m_error_abs = [m_error_abs[0], m_error_abs[1], m_error_abs[2],
                   m_error_abs[4], m_error_abs[5], m_error_abs[6],]

    ax = plt.subplot(G[2, :])
    plot1 = ax.boxplot(m_error_abs,
                       positions=tem_x,
                       medianprops={'color': 'orange'},
                       whis=(0, 100),
                       boxprops=dict(facecolor=c1[1],  # 箱体填充颜色
                                     color='k',  # 箱体边界颜色
                                     linewidth=0.5,  # 箱体边界粗细
                                     linestyle='-'),
                       capprops=dict(color=c1[5],  # 箱体边界颜色
                                     linewidth=1.,  # 箱体边界粗细
                                     linestyle='-'),
                       patch_artist=True,
                       showfliers=False,
                       widths=.28,
                       notch=False)
    for patch in plot1['boxes']:
        patch.set_facecolor(c1[1])
    [patch.set(alpha=0.6) for patch in plot1['boxes']]  #####设置颜色的透明度

    plt.ylabel(r'$Error\ percentage\ (\%)$', fontsize=lsize)
    plt.tick_params(labelsize=lsize)
    plt.title('d', x=-0.064, y=1.1, font=font_l, horizontalalignment='left', verticalalignment='top')
    plt.grid(linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
    plt.xticks(tem_x, tem_x_label, fontsize=lsize)
    plt.ylim([0, 2.5])

    plt.tight_layout()
    plt.show()





    time_end = Time.time()
    print('cost %f s' % (time_end - time_start))


