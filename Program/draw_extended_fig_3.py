import netCDF4 as nc
import numpy as np
import scipy.io as scio
import warnings
import time as Time
import matplotlib
import matplotlib.pyplot as plt
import scipy.ndimage
import matplotlib.gridspec as gridspec
import cmaps
import scienceplots
import matplotlib.colors as mcolors


if __name__=='__main__':
    time_begin=Time.time()
    warnings.filterwarnings('ignore')
    matplotlib.rc("font", family='Verdana', weight="normal")
    matplotlib.rc('text', usetex=True)
    matplotlib.style.use(['science', 'nature', 'no-latex'])


    norm = mcolors.TwoSlopeNorm(vmin=-12., vcenter=0, vmax=5.)
    m_cmap = cmaps.MPL_bwr
    lsize = 18
    fig = plt.figure(figsize=(16.5, 8.5), dpi=200)  # 设置大小和分辨率
    font = {'family': 'Verdana', 'weight': 'bold', 'size': lsize + 2, }
    font_l = {'family': 'Verdana', 'weight': 'normal', 'size': lsize, }
    font_l1 = {'family': 'Verdana', 'weight': 'bold', 'size': lsize, }

    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_5/anomaly.mat')
    anomaly_pco2 = tem_data['anomaly_pco2']
    anomaly_sst = tem_data['anomaly_sst']
    anomaly_sss = tem_data['anomaly_sss']
    anomaly_dic = tem_data['anomaly_dic']
    anomaly_ta = tem_data['anomaly_ta']


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    G = gridspec.GridSpec(5, 4, figure=fig)
    xx, yy = np.mgrid[-10:10:201j, -10:10:201j]
    rmax = 1
    m_lim = [-7.2, 7.2]
    ll = 20
    rr = 181
    al=1
    lsize1 = lsize + 1.5
    m_lin = np.linspace(-12, 5, 18)

    ax = plt.subplot(G[0:4, 0:2])
    con1 = plt.contourf(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_pco2[ll:rr, ll:rr],
                 m_lin, cmap=m_cmap, alpha=al, extend='both', norm=norm)
    con_line = plt.contour(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_pco2[ll:rr, ll:rr],
                        [-7,-6,-5,-4],
                        colors='w', linewidth=3, alpha=1)

    plt.clabel(con_line, fmt='%.1f', fontsize=lsize1 + 2, inline=0.01, inline_spacing=0.4, colors='k')
    ax.tick_params(labelsize=lsize)
    plt.hlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.vlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.xlim(m_lim)
    plt.ylim(m_lim)
    plt.ylabel(r'Distance along track ($R_{\rm max}$)', fontsize=lsize)
    ax.yaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.quiver(0, 0, 0, 1, color='w', angles='xy', scale=3., width=0.02, headwidth=3.5,headlength=4.9)
    plt.title(r'a  Surface $\bf{pCO_2^{sea}}$ anomaly', x=0.02, y=0.98,  font=font, horizontalalignment='left',
              verticalalignment='top')

    ax = plt.subplot(G[0:2, 2:3])
    con = plt.contourf(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_sst[ll:rr, ll:rr],
                 m_lin, cmap=m_cmap, alpha=al, extend='both', norm=norm)
    con_line = plt.contour(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_sst[ll:rr, ll:rr],
                        [-11,-10,-9,-8,-7],
                        colors='w', linewidth=1, alpha=0.75)
    plt.clabel(con_line, fmt='%.1f', fontsize=lsize1, inline=0.01, inline_spacing=0.4, colors='k')
    ax.tick_params(labelsize=lsize)
    plt.xlim(m_lim)
    plt.ylim(m_lim)
    plt.xticks([-5, -2.5, 0, 2.5, 5], [], fontsize=lsize)
    plt.hlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.vlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.text(0.03, 0.98, r'b SST contribution', font=font, rotation='horizontal',
             horizontalalignment='left', verticalalignment='top',transform=ax.transAxes)


    ax = plt.subplot(G[0:2, 3:4])
    con = plt.contourf(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_sss[ll:rr, ll:rr],
                 m_lin, cmap=m_cmap, alpha=al, extend='both', norm=norm)
    con_line = plt.contour(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_sss[ll:rr, ll:rr],
                        [-1, 0, 1],
                        colors='w', linewidth=1, alpha=0.75)
    plt.clabel(con_line, fmt='%.1f', fontsize=lsize1, inline=0.01, inline_spacing=0.4, colors='k')
    ax.tick_params(labelsize=lsize)
    plt.xlim(m_lim)
    plt.ylim(m_lim)
    plt.xticks([-5, -2.5, 0, 2.5, 5], [], fontsize=lsize)
    plt.yticks([-5, -2.5, 0, 2.5, 5], [], fontsize=lsize)
    plt.hlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.vlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.text(0.03, 0.98, r'c SSS contribution', font=font, rotation='horizontal',
             horizontalalignment='left', verticalalignment='top',transform=ax.transAxes)

    ax = plt.subplot(G[2:4, 2:3])
    con = plt.contourf(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_dic[ll:rr, ll:rr],
                 m_lin, cmap=m_cmap, alpha=al, extend='both', norm = norm)
    con_line = plt.contour(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_dic[ll:rr, ll:rr],
                        [1,2,3,4,5],
                        colors='w', linewidth=1, alpha=0.75)
    plt.clabel(con_line, fmt='%.1f', fontsize=lsize1, inline=0.01, inline_spacing=0.4, colors='k')
    ax.tick_params(labelsize=lsize)
    plt.xlim(m_lim)
    plt.ylim(m_lim)

    ax.yaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.hlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.vlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.text(0.03, 0.98, r'd DIC contribution', font=font, rotation='horizontal',
             horizontalalignment='left', verticalalignment='top',transform=ax.transAxes)


    ax = plt.subplot(G[2:4, 3:4])
    con = plt.contourf(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_ta[ll:rr, ll:rr],
                 m_lin, cmap=m_cmap, alpha=al, extend='both', norm=norm)
    con_line = plt.contour(xx[ll:rr, ll:rr] / rmax, yy[ll:rr, ll:rr] / rmax, anomaly_ta[ll:rr, ll:rr],
                        [0,1,2,3,4],
                        colors='w', linewidth=1, alpha=0.75)
    plt.clabel(con_line, fmt='%.1f', fontsize=lsize1, inline=0.01, inline_spacing=0.4, colors='k')
    ax.tick_params(labelsize=lsize)
    plt.xlim(m_lim)
    plt.ylim(m_lim)

    ax.yaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.xaxis.set_major_locator(plt.MultipleLocator(2.5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    plt.yticks([-5, -2.5, 0, 2.5, 5], [], fontsize=lsize)
    plt.hlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.vlines(0, -20, 20, colors='k', linestyles='--', linewidth=1.0, alpha=0.9)
    plt.text(0.03, 0.98, r'e TA contribution', font=font, rotation='horizontal',
             horizontalalignment='left', verticalalignment='top',transform=ax.transAxes)

    fig.subplots_adjust(top=0.9)
    cbar_ax = fig.add_axes([0.1, 0.069, 0.8, 0.04])
    cbar = fig.colorbar(con, cax=cbar_ax, orientation='horizontal', fraction=.057, pad=0.03)
    cbar.set_ticks(np.linspace(-12,5,18))
    cbar.ax.tick_params(labelsize=lsize-1)
    cbar.ax.set_xlabel(r'$pCO_2^{sea}$ anomaly ($\mu$atm)', fontsize=lsize-1)

    ax = plt.subplot(G[4, 0:2])
    plt.xlim(m_lim)
    plt.text(0, 0.9, r'Distance across track ($R_{\rm max}$)', fontsize=lsize - 1, rotation='horizontal',
             horizontalalignment='center', verticalalignment='center')
    plt.axis('off')
    ax = plt.subplot(G[4, 2:3])
    plt.xlim(m_lim)
    plt.text(0, 0.9, r'Distance across track ($R_{\rm max}$)', fontsize=lsize - 1, rotation='horizontal',
             horizontalalignment='center', verticalalignment='center')
    plt.axis('off')
    ax = plt.subplot(G[4, 3:4])
    plt.xlim(m_lim)
    plt.text(0, 0.9, r'Distance across track ($R_{\rm max}$)', fontsize=lsize - 1, rotation='horizontal',
             horizontalalignment='center', verticalalignment='center')
    plt.axis('off')

    plt.tight_layout()
    plt.show()


    time_end = Time.time()
    print('cost %f s' % (time_end-time_begin))




