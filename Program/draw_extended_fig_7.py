import warnings
import time as Time
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize, LinearSegmentedColormap
import numpy as np
import scipy.io as scio
from scipy.optimize import curve_fit
import uncertainties as unc
from scipy.stats import norm
import scienceplots
import seaborn as sn



if __name__ == '__main__':
    time_begin = Time.time()
    warnings.filterwarnings('ignore')
    matplotlib.rc("font", family='Verdana', weight="normal")
    matplotlib.style.use(['science', 'nature', 'no-latex'])



    cmap_modify_add  = [sn.light_palette('#61BDCB', n_colors=10)[-1]]
    cmap_modify_add.append(sn.light_palette('#50BBD5', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#3EAFD9', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#478BBF', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#3F6EAA', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#395FA0', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#34569B', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#455B9D', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#6E78AF', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#A5A6CC', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#DDDDED', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#FADEDC', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#EEAAAB', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E78183', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E1605F', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E05248', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E36845', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E97E44', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#EE9847', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#F2B74A', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#F2D754', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#E5E051', n_colors=10)[-1])
    cmap_modify_all = LinearSegmentedColormap.from_list('mycmap', cmap_modify_add, N=30)



    xx, yy = np.mgrid[0:359.75:1440j, -90 + 0.125:90 - 0.125:720j]
    xx[xx < 30] = xx[xx < 30] + 360
    xx = np.vstack((xx[120:, :], xx[:120, :]))


    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_11/Mix_pco2_trend_mean.mat')
    Mix_pco2_trend_mean_126 = tem_data['Mix_pco2_trend_mean_126'] * 10. # change the unit from uatm/yr to uatm/dec
    Mix_pco2_trend_mean_585 = tem_data['Mix_pco2_trend_mean_585'] * 10. # change the unit from uatm/yr to uatm/dec

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    lsize = 20
    m_al=0.2
    font = {'family': 'Verdana', 'weight': 'bold', 'size': lsize + 2.5, }
    font_l = {'family': 'Verdana', 'weight': 'normal', 'size': lsize - 1, }
    lw = 2.
    al = 1.

    fig = plt.figure(figsize=(11, 9), dpi=300)  # 设置大小和分辨率
    G = gridspec.GridSpec(7, 4, figure=fig)
    ax = plt.subplot(G[0:3, :])
    lsize = 15
    lllon = 30
    lllat = -60
    urlon = 390.01
    urlat = 60.01

    m = Basemap(projection='cyl', resolution='l', llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat)
    m.drawcoastlines(color='gray', linewidth=0.4, )
    m.drawmapboundary(fill_color='white')

    m.fillcontinents(color='gray', alpha=1, lake_color='gray')
    m.drawmeridians(np.arange(lllon, urlon, 30), labels=[0, 0, 0, 1], fontsize=lsize, linewidth=0.8,
                    color='gray')  # 经线
    m.drawparallels(np.arange(lllat, urlat, 20), labels=[1, 0, 0, 0], fontsize=lsize, linewidth=0.8,
                    color='gray')  # 纬线
    tem_data = Mix_pco2_trend_mean_126 * 1.
    plt.contourf(xx, yy, tem_data, np.linspace(-10, 10, 61), colors='white', extend='both', alpha=0.9)
    plt.gca().set_facecolor("#C0C0C0")
    con = plt.contourf(xx, yy, tem_data, np.linspace(-10, 10, 61),
                       cmap=cmap_modify_all, extend='both', alpha=0.9, )


    m_lw=2.
    m_ls='--'
    m_lc='k'
    m_la=0.6
    l = 100
    r = 180
    u = 40
    d = 10
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 180
    r = 255
    u = 60
    d = 10
    m.plot([r, r], [30, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, 263], [30, 18], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([263, 275], [18, 12.2], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([275, 277], [12.2, 8.6], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, 260], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 5, 30, r'NEP', font=font_l, horizontalalignment='right', verticalalignment='center', )
    m.plot([180, 180], [10, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(100 + 5, 30, r'NWP', font=font_l, horizontalalignment='left', verticalalignment='center', )
    m.plot([100, 360], [40, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 260
    r = 359
    u = 40
    d = 10
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 5, 30, r'NA', font=font_l, horizontalalignment='right', verticalalignment='center', )

    l = 31
    r = 135
    u = -10
    d = -40
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 3, u - 3, r'SI', font=font_l, horizontalalignment='right', verticalalignment='top', )

    l = 135
    r = 240
    u = -10
    d = -40
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 3, u - 3, r'SP', font=font_l, horizontalalignment='right', verticalalignment='top', )
    m.plot([30, 240], [-40, -40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 33
    r = 100
    u = 30
    d = 10
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(47, 21, r'NI', font=font_l, horizontalalignment='right', verticalalignment='center', )
    plt.title('a', x=-0.11, y=1,  font=font, horizontalalignment='left',
              verticalalignment='top')



    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    ax = plt.subplot(G[3:6, :])
    lsize = 17
    lllon = 30
    lllat = -60
    urlon = 390.01
    urlat = 60.01

    m = Basemap(projection='cyl', resolution='l', llcrnrlon=lllon, llcrnrlat=lllat, urcrnrlon=urlon, urcrnrlat=urlat)
    m.drawcoastlines(color='gray', linewidth=0.4, )
    m.drawmapboundary(fill_color='white')

    m.fillcontinents(color='gray', alpha=1, lake_color='gray')
    m.drawmeridians(np.arange(lllon, urlon, 30), labels=[0, 0, 0, 1], fontsize=lsize, linewidth=0.8,
                    color='gray')  # 经线
    m.drawparallels(np.arange(lllat, urlat, 20), labels=[1, 0, 0, 0], fontsize=lsize, linewidth=0.8,
                    color='gray')  # 纬线
    tem_data = Mix_pco2_trend_mean_585 * 1.
    plt.contourf(xx, yy, tem_data, np.linspace(-10, 10, 61), colors='white', extend='both', alpha=0.9)
    plt.gca().set_facecolor("#C0C0C0")
    con = plt.contourf(xx, yy, tem_data, np.linspace(-10, 10, 61),
                       cmap=cmap_modify_all, extend='both', alpha=0.9, )
    lw = 2.9
    m_c = 'k'
    m_al = 0.6

    m_lw=2.
    m_ls='--'
    m_lc='k'
    m_la=0.6
    l = 100
    r = 180
    u = 40
    d = 10
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 180
    r = 255
    u = 60
    d = 10
    m.plot([r, r], [30, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, 263], [30, 18], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([263, 275], [18, 12.2], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([275, 277], [12.2, 8.6], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, 260], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 5, 30, r'NEP', font=font_l, horizontalalignment='right', verticalalignment='center', )
    m.plot([180, 180], [10, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(100 + 5, 30, r'NWP', font=font_l, horizontalalignment='left', verticalalignment='center', )
    m.plot([100, 360], [40, 40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 260
    r = 359
    u = 40
    d = 10
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 5, 30, r'NA', font=font_l, horizontalalignment='right', verticalalignment='center', )

    l = 31
    r = 135
    u = -10
    d = -40
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 3, u - 3, r'SI', font=font_l, horizontalalignment='right', verticalalignment='top', )

    l = 135
    r = 240
    u = -10
    d = -40
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(r - 3, u - 3, r'SP', font=font_l, horizontalalignment='right', verticalalignment='top', )
    m.plot([30, 240], [-40, -40], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)

    l = 33
    r = 100
    u = 30
    d = 10
    m.plot([l, l], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([r, r], [d, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [u, u], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    m.plot([l, r], [d, d], lw=m_lw, c=m_lc, ls=m_ls, alpha=m_la)
    ax.text(47, 21, r'NI', font=font_l, horizontalalignment='right', verticalalignment='center', )
    plt.title('b', x=-0.11, y=1,  font=font, horizontalalignment='left',
              verticalalignment='top')



    fig.subplots_adjust(top=0.9)
    cbar_ax = fig.add_axes([0.12, 0.075, 0.82, 0.035])
    cbar = fig.colorbar(con, cax=cbar_ax, orientation='horizontal', fraction=.057, pad=0.05)
    cbar.set_ticks(np.linspace(-10, 10, 11))
    cbar.ax.tick_params(labelsize=lsize)
    cbar.ax.set_xlabel(r'Trend ($\mu$atm $\rm dec^{-1}$)', fontsize=lsize - 1)
    plt.tight_layout()
    plt.show()



    time_end = Time.time()
    print('cost %f s' % (time_end-time_begin))

