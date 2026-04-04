import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage
from mpl_toolkits.basemap import Basemap
import matplotlib
from matplotlib.colors import Normalize,LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import scipy.ndimage
import scipy.io as scio
import warnings
import time as Time
import seaborn as sn
import scienceplots


if __name__=='__main__':
    time_begin=Time.time()
    warnings.filterwarnings('ignore')
    matplotlib.style.use(['science', 'nature', 'no-latex'])

    c2 = ['#01847F', '#703274', ]

    cmap_modify_add = [sn.light_palette('#f1f1f1', n_colors=10)[-1]]
    cmap_modify_add.append(sn.light_palette('#f1f1f1', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#DDDDED', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#b5b7d7', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#8489ba', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#5463a2', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#37549a', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#255097', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#274892', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#463789', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#752875', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#a6214e', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#db2d31', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#de4134', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#e15936', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#e6783c', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#ee9c43', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#f5ca4c', n_colors=10)[-1])
    cmap_modify_add.append(sn.light_palette('#e6e153', n_colors=10)[-1])
    cmap_modify_all = LinearSegmentedColormap.from_list('mycmap', cmap_modify_add, N=20)



    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

    ew = 2
    lw = 4
    lsize = 19.5
    m_al=0.9
    m_capsize = 4.5
    m_capthick = 2.5
    font = {'family': 'Verdana', 'weight': 'bold', 'size': lsize + 4, }
    lw = 3.


    fig=plt.figure(figsize=(15., 23), dpi=200)
    G = gridspec.GridSpec(29, 4, figure=fig)
    ax = plt.subplot(G[0:8,0:2])
    tem_data = scio.loadmat('./intermediate_file/Fig_1/method_LIDR.mat')
    x_range = tem_data['x_range'][0]
    y_range = tem_data['y_range'][0]
    bins_DIC = tem_data['bins_DIC']

    lin_DIC=np.linspace(y_range[0], y_range[1], 7)
    num = 120
    lin_space = (x_range[1] - x_range[0]) / num / 2
    lin = np.linspace(0, 2, 11)

    x, y = np.mgrid[x_range[0] :x_range[1] :(num+1) * 1j, y_range[0]:y_range[1]:(num+1) * 1j]
    con = ax.pcolormesh(x, y, np.log10(bins_DIC),
                        vmax=2, vmin=0,
                        cmap=cmap_modify_all, antialiased=True
                        )
    ax.plot(x_range, y_range, linewidth=1.2, c='k')

    plt.grid(linestyle='--',linewidth=0.6,c='gray',alpha=0.8)
    ax.tick_params(labelsize=lsize)
    ax.text(2050, 1825, r'$\ R^2 \ \ \ \ \ =\ '+'%.3f' % 0.943+'$', horizontalalignment='left',fontsize=lsize-2)
    ax.text(2050, 1780, r'$\ \ n \ \ \ \ \ \ =\ '+'%d' % 11058+'$', horizontalalignment='left', fontsize=lsize-2)
    ax.text(2050, 1735, r'$RMSE\ =\ '+ '%.1f' % 17.6 +'$', horizontalalignment='left', fontsize=lsize-2)
    plt.xticks([1700,1900,2100,2300],[1700,1900,2100,2300],fontsize=lsize)
    plt.yticks([1700, 1900, 2100, 2300], [1700, 1900, 2100, 2300], fontsize=lsize)
    plt.xlim([1720, 2280])
    plt.ylim([1720, 2280])
    ax.set_ylabel(r'${\rm DIC_{Esti}}\ ({\mu} \rm mol\ kg^{-1})$', fontsize=lsize)
    plt.title('a', x=-0.12, y=0.99, font=font, horizontalalignment='left', verticalalignment='top')
    fig.subplots_adjust(top=0.9)
    cbar_ax = fig.add_axes([0.10, 0.96, 0.25, 0.017])
    cbar = fig.colorbar(con, cax=cbar_ax, orientation='horizontal', fraction=.057, pad=0.05)
    cbar.set_ticks([0, 0.5, 1, 1.5, 2])
    cbar.ax.tick_params(labelsize=lsize - 2)
    cbar.set_label(r'$Log_{10}\ frequency$', fontsize=lsize - 2)


    ax = plt.subplot(G[0:8,2:4])
    tem_data = scio.loadmat('./intermediate_file/Fig_1/method_LIAR.mat')
    x_range = tem_data['x_range'][0]
    y_range = tem_data['y_range'][0]
    bins_TA = tem_data['bins_TA']
    lin_TA=np.linspace(y_range[0],y_range[1],7)
    num = 120
    lin_space = (x_range[1] - x_range[0]) / num / 2
    lin = np.linspace(0, 2, 11)
    x, y = np.mgrid[x_range[0] :x_range[1] :(num+1) * 1j, y_range[0]:y_range[1]:(num+1) * 1j]
    con = ax.pcolormesh(x, y, np.log10(bins_TA),
                        vmax=2, vmin=0,
                        cmap=cmap_modify_all, antialiased=True
                        )
    ax.plot(x_range, y_range, linewidth=1.2, c='k')

    plt.grid(linestyle='--',linewidth=0.6,c='gray',alpha=0.8)
    ax.tick_params(labelsize=lsize)
    ax.text(2350, 2125, r'$\ R^2\ \ \ \ \ =\ '+'%.3f' % 0.978+'$', horizontalalignment='left',fontsize=lsize-2)
    ax.text(2350, 2080, r'$\ \ n\ \ \ \ \ \ =\ '+'%d' % 9347+'$', horizontalalignment='left', fontsize=lsize-2)
    ax.text(2350, 2035, r'$RMSE\ =\ '+ '%.1f' % 8.7 +'$', horizontalalignment='left', fontsize=lsize-2)
    plt.xticks([2000, 2200, 2400, 2600], [2000, 2200, 2400, 2600],fontsize=lsize)
    plt.yticks([2000, 2200, 2400, 2600], [2000, 2200, 2400, 2600], fontsize=lsize)
    plt.xlim([2020, 2580])
    plt.ylim([2020, 2580])
    ax.set_ylabel(r'${\rm TA_{Esti}}\ ({\mu} \rm mol\ kg^{-1})$', fontsize=lsize)
    plt.title('b', x=-0.12, y=0.99, font=font, horizontalalignment='left',
                 verticalalignment='top')


    ax = plt.subplot(G[8, :])
    plt.text(0.225, 1, r'${\rm DIC_{Meas}}\ ({\mu} \rm mol\ kg^{-1})$', fontsize=lsize, horizontalalignment='center')
    plt.text(0.775, 1, r'${\rm TA_{Meas}}\ ({\mu} \rm mol\ kg^{-1})$', fontsize=lsize, horizontalalignment='center')
    plt.axis('off')




    #==================================================================================================================#
    ax = plt.subplot(G[9:15, :])
    tem_data = scio.loadmat('./intermediate_file/Fig_1/evo_dpco2.mat')
    dpco2_obs = np.squeeze(tem_data['dpco2_obs'])
    dpco2_rec = np.squeeze(tem_data['dpco2_rec'])
    dpco2_obs_std = np.squeeze(tem_data['dpco2_obs_std'])
    dpco2_rec_std = np.squeeze(tem_data['dpco2_rec_std'])

    l1, = plt.plot(np.linspace(-20, 59, 80), dpco2_obs, linewidth=lw, c=c2[0], linestyle='-', alpha=m_al, )
    l2, = plt.plot(np.linspace(-20, 59, 80), dpco2_rec, linewidth=lw, c=c2[1], linestyle='-', alpha=m_al, )
    tem_error = dpco2_obs_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], dpco2_obs[::2], yerr=tem_error[::2], fmt='',
                 ecolor=c2[0],linestyle='none', alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    tem_error = dpco2_rec_std
    plt.errorbar(np.linspace(-20, 59, 80)[::2], dpco2_rec[::2], yerr=tem_error[::2], fmt='',
                 ecolor=c2[1],linestyle='none', alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    plt.xticks(np.linspace(-20, 60, 9), [])
    plt.xlim([-25,65])
    plt.ylim([-14,4])
    plt.yticks([-10, -5, 0])
    ax.tick_params(labelsize=lsize)
    plt.ylabel(r'$DpCO_2\ (\mu atm)$', fontsize=lsize)
    plt.fill_betweenx([-20, 10], -10, -3, alpha=0.1, color='gray')
    plt.fill_betweenx([-20, 10], -3, 3, alpha=0.1, color='gray')
    plt.fill_betweenx([-20, 10], 3, 65, alpha=0.2, color='gray')
    plt.vlines(0, -14, 4, colors='k', alpha=0.5, linewidth = 2.0, linestyles='--')
    tem_data_1 = dpco2_obs * 1.
    tem_data_2 = dpco2_rec * 1.
    tem_r = np.corrcoef(tem_data_1, tem_data_2)[0,1]
    tem_bias = - np.nanmean(tem_data_1 - tem_data_2)
    ax.text(46, -12,r'$\ \ r\ \ \ = $' + '%.02f' % tem_r,  fontsize=lsize-2, c='k', alpha=0.8)
    ax.text(46, -13, r'$bias =' + '%.02f' % tem_bias+' $', fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(-23.5, -12, r'$DpCO_2 = pCO_2^{sea} - pCO_2^{air}$', fontsize=lsize - 2, c='k', alpha=0.8, horizontalalignment='left',)
    ax.text(31, 1, r'$Recovery$', fontsize=lsize - 2, c='k', alpha=1, horizontalalignment='center', verticalalignment='center')
    plt.title('c', x=-0.06, y=0.99, font=font, horizontalalignment='left', verticalalignment='top')
    plt.legend([l1, l2], [r'$DpCO_{2}^{obs}$', r'$DpCO_{2}^{rec}$'], fontsize=lsize - 2, loc=1, framealpha=0.6,
               edgecolor='gray', frameon=False, ncol=2)
    ax.text(-23.5, 2, r'$37\ TCs$', fontsize=lsize, c='k', alpha=0.8, horizontalalignment='left',)




    #==================================================================================================================#
    ax = plt.subplot(G[15:21, :])
    tem_data = scio.loadmat('./intermediate_file/Fig_1/evo_wind.mat')
    wind_obs = np.squeeze(tem_data['wind_obs'])
    wind_rec = np.squeeze(tem_data['wind_rec'])
    wind_obs_std = np.squeeze(tem_data['wind_obs_std'])
    wind_rec_std = np.squeeze(tem_data['wind_rec_std'])

    plt.plot(np.linspace(-20, 59, 80), wind_rec, linewidth=lw, linestyle='-', c=c2[1], alpha=m_al)
    plt.plot(np.linspace(-20, 59, 80), wind_obs, linewidth=lw, linestyle='-', c=c2[0], alpha=m_al )
    tem_error = wind_obs_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], wind_obs[::2], yerr=tem_error[::2],
                 fmt='',linestyle='none',
                 ecolor=c2[0], alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    tem_error = wind_rec_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], wind_rec[::2], yerr=tem_error[::2], fmt='',
                 ecolor=c2[1],linestyle='none',alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.xaxis.set_major_locator(plt.MultipleLocator(10))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))
    ax.tick_params(labelsize=lsize)
    plt.xticks(np.linspace(-20, 60, 9), [])
    plt.fill_betweenx([0, 20], -10, -3, alpha=0.1, color='gray')
    plt.fill_betweenx([0, 20], -3, 3, alpha=0.2, color='gray')
    plt.fill_betweenx([0, 20],  3, 65, alpha=0.1, color='gray')
    plt.vlines(0, 4, 16, colors='k', alpha=0.5, linewidth = 2.0, linestyles='--')
    plt.xlim([-25, 65])
    plt.ylim([4,16])

    tem_data_1 = wind_obs * 1.
    tem_data_2 = wind_rec * 1.
    tem_r = np.corrcoef(tem_data_1, tem_data_2)[0,1]
    tem_bias = - np.nanmean(tem_data_1 - tem_data_2)
    ax.text(46, 11, r'$\ \ r\ \ \ =' + '%.02f' % tem_r + ' $', fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(46, 10, r'$bias =' + '%.02f' % tem_bias + ' $', fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(0, 5, r'$Forcing$', fontsize=lsize - 2, c='k', alpha=1, horizontalalignment='center', verticalalignment='center')
    plt.ylabel(r'$v_{10}\ (m\ s^{-1})$', fontsize=lsize)
    plt.title('d', x=-0.06, y=0.99, font=font, horizontalalignment='left', verticalalignment='top')
    plt.legend([l1, l2], [r'$v_{10}^{obs}$', r'$v_{10}^{rec}$'], fontsize=lsize - 2, loc=1, framealpha=0.6,
               edgecolor='gray', frameon=False, ncol=2)





    #==================================================================================================================#
    ax = plt.subplot(G[21:29, :])
    tem_data = scio.loadmat('./intermediate_file/Fig_1/evo_fco2_positive.mat')
    fco2_obs = np.squeeze(tem_data['fco2_obs'])
    fco2_rec = np.squeeze(tem_data['fco2_rec'])
    fco2_obs_std = np.squeeze(tem_data['fco2_obs_std'])
    fco2_rec_std = np.squeeze(tem_data['fco2_rec_std'])


    l1, = plt.plot(np.linspace(-20, 59, 80), fco2_obs, linewidth=lw, c=c2[0], linestyle='-', alpha=m_al, )
    l2, = plt.plot(np.linspace(-20, 59, 80), fco2_rec, linewidth=lw, c=c2[1], linestyle='-', alpha=m_al, )
    tem_error = fco2_obs_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], fco2_obs[::2], yerr=tem_error[::2],
                 fmt='', ecolor=c2[0], linestyle='none', alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    tem_error = fco2_rec_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], fco2_rec[::2], yerr=tem_error[::2],
                 fmt='', linestyle='none',
                 ecolor=c2[1], alpha=0.8,
                 elinewidth=lw-0.5, capsize=m_capsize, capthick=m_capthick)
    tem_data_1 = fco2_obs * 1.
    tem_data_2 = fco2_rec * 1.
    tem_r = np.corrcoef(tem_data_1, tem_data_2)[0,1]
    tem_bias = - np.nanmean(tem_data_1 - tem_data_2)
    ax.text(-23, 11, r'$\ \ r\ \ \ = %.02f$' % tem_r, fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(-23, 9, r'$bias =' + '%.02f' % tem_bias + ' $', fontsize=lsize - 2, c='k', alpha=0.8)
    #------------------------------------------------------------------------------------------------------------------#
    tem_data = scio.loadmat('./intermediate_file/Fig_1/evo_fco2_negative.mat')
    fco2_obs = np.squeeze(tem_data['fco2_obs'])
    fco2_rec = np.squeeze(tem_data['fco2_rec'])
    fco2_obs_std = np.squeeze(tem_data['fco2_obs_std'])
    fco2_rec_std = np.squeeze(tem_data['fco2_rec_std'])


    l1, = plt.plot(np.linspace(-20, 59, 80), fco2_obs, linewidth=lw, c=c2[0], linestyle='-', alpha=m_al, )
    l2, = plt.plot(np.linspace(-20, 59, 80), fco2_rec, linewidth=lw, c=c2[1], linestyle='-', alpha=m_al, )
    tem_error = fco2_obs_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], fco2_obs[::2], yerr=tem_error[::2],
                 fmt='',
                 ecolor=c2[0], linestyle='none', alpha=0.8,
                 elinewidth=lw - 0.5, capsize=m_capsize, capthick=m_capthick)
    tem_error = fco2_rec_std * 1.
    plt.errorbar(np.linspace(-20, 59, 80)[::2], fco2_rec[::2], yerr=tem_error[::2],
                 fmt='', linestyle='none',
                 ecolor=c2[1], alpha=0.8,
                 elinewidth=lw - 0.5, capsize=m_capsize, capthick=m_capthick)

    plt.xlabel(r'$Time\ (d)$', fontsize=lsize)
    plt.ylabel(r'$CO_2\ flux\ (mmol\ m^{-2}\ d^{-1})$', fontsize=lsize)
    ax.tick_params(labelsize=lsize)
    plt.xlim([-25, 65])
    plt.legend([l1, l2], [r'$F_{CO_2}^{obs}$', r'$F_{CO_2}^{rec}$'], fontsize=lsize-2, loc=4, framealpha=0.6, edgecolor='gray', frameon=False, ncol=2)

    tem_data_1 = fco2_obs * 1.
    tem_data_2 = fco2_rec * 1.
    tem_r = np.corrcoef(tem_data_1, tem_data_2)[0,1]
    tem_bias = - np.nanmean(tem_data_1 - tem_data_2)
    ax.text(-23, -10, r'$\ \ r\ \ \ = %.02f$' % tem_r, fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(-23, -12, r'$bias =' + '%.02f' % tem_bias + ' $', fontsize=lsize - 2, c='k', alpha=0.8)
    ax.text(-6.5, 10, r'$Prestorm$', fontsize=lsize - 2, c='k', alpha=1,
            horizontalalignment='center', verticalalignment='center')

    plt.vlines(0, -20, 20, colors='k', alpha=0.5, linewidth = 2.0, linestyles='--')
    plt.hlines(0, -25, 65, colors='k', alpha=0.5, linewidth = 2.0, linestyles='-')
    plt.fill_betweenx([-20, 20], -10, -3, alpha=0.2, color='gray')
    plt.fill_betweenx([-20, 20], -3, 3, alpha=0.1, color='gray')
    plt.fill_betweenx([-20, 20],  3, 65, alpha=0.1, color='gray')
    ax.text(-23.5, 16, r'$19\ TCs$', fontsize=lsize, c='k', alpha=0.8, horizontalalignment='left', )
    ax.text(-23.5, -16, r'$18\ TCs$', fontsize=lsize, c='k', alpha=0.8, horizontalalignment='left',  verticalalignment='top')

    ax.text(10, 14, r'$DpCO_2\ >\ 0$', fontsize=lsize, c='k', alpha=0.8, horizontalalignment='left', )
    ax.text(10, -14, r'$DpCO_2\ <\ 0$', fontsize=lsize, c='k', alpha=0.8, horizontalalignment='left',  verticalalignment='top')

    plt.ylim([-19, 19])
    plt.title('e', x=-0.06, y=0.99, font=font, horizontalalignment='left', verticalalignment='top')



    #==================================================================================================================#
    ax = ax.inset_axes([0.62,0.61,0.36,0.37])

    m = Basemap(ax=ax, projection='cyl', resolution='l', llcrnrlon=30, llcrnrlat=-75, urcrnrlon=390, urcrnrlat=75.001)
    m.drawmapboundary(fill_color='white')
    m.fillcontinents(color='gray', alpha=0.6, lake_color='gray')
    m.drawmeridians(np.arange(30, 390.001, 60), labels=[0, 0, 0, 1], fontsize=lsize-3, linewidth=0.8,
                    color='gray')  # 经线
    m.drawparallels(np.arange(-60, 60.001, 30), labels=[1, 0, 0, 0], fontsize=lsize-3, linewidth=0.8,
                    color='gray')  # 纬线
    ls=45
    ts=20
    ax.scatter(144.58,32.28, color='r',marker='o',s=ls)
    ax.scatter(90, 15, color='r', marker='o', s=ls)
    ax.scatter(360-157.98, 22.67, color='r', marker='o', s=ls)
    ax.scatter(360-64.2, 31.5, color='r', marker='o', s=ls)
    ax.scatter(360-12.67, 68, color='r', marker='o', s=ls)
    ax.text(144.58,32.28 - ts, 'KEO', horizontalalignment='center', verticalalignment='center', fontsize=lsize - 4.)
    ax.text(360-157.98, 22.67 - ts, 'WHOTS', horizontalalignment='center', verticalalignment='center', fontsize=lsize - 4.)
    ax.text(360-64.2, 31.5 - ts, 'BTM', horizontalalignment='center', verticalalignment='center', fontsize=lsize - 4.)
    ax.text(360-12.67, 68 -ts, 'ICELAND', horizontalalignment='center', verticalalignment='center', fontsize=lsize - 4.)
    ax.text(90, 15 - ts, 'BOBOA', horizontalalignment='center', verticalalignment='center', fontsize=lsize - 4.)

    ax.tick_params(labelsize=lsize-2)


    plt.tight_layout()
    G.tight_layout(fig, pad=0.4, h_pad=None, w_pad=None, rect=None)
    plt.show()


    time_end = Time.time()
    print('cost %f s' % (time_end-time_begin))


