import netCDF4 as nc
import numpy as np
import scipy.io as scio
import warnings
import time as Time
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import scipy.ndimage
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm
from matplotlib.colors import Normalize,LinearSegmentedColormap
from sklearn.linear_model import LinearRegression
from scipy.interpolate import interp1d
import os
import cmaps
import scienceplots
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit


def m_fit1d(x, y, xi, e, error_on):
    f = lambda x, *p: np.polyval(p, x)  # 定义待拟合函数
    # 依次传入拟合函数、拟合数据，参数初值，及误差
    if error_on:
        popt, cov = curve_fit(f, x, y, [1, 1], sigma=e, absolute_sigma=True)
    else:
        popt, cov = curve_fit(f, x, y, [1, 1], sigma=None, absolute_sigma=False)
    ps = np.random.multivariate_normal(popt, cov, 20000)
    ysample = np.asarray([f(xi, *pi) for pi in ps])
    lower = np.percentile(ysample, 2.5, axis=0)
    upper = np.percentile(ysample, 97.5, axis=0)
    m_line = np.poly1d(popt)(xi)
    return lower, upper, m_line, popt, cov

if __name__=='__main__':

    time_begin=Time.time()
    warnings.filterwarnings('ignore')
    matplotlib.rc("font", family='Verdana', weight="normal")
    matplotlib.rc('text', usetex=True)
    matplotlib.style.use(['science', 'nature', 'no-latex'])



    c2 = ['#01847F', '#703274', '#1762a1', ]


    #==================================================================================================================#

    al = 1.
    m_al = 0.2
    lsize = 16
    font_l = {'family': 'Verdana', 'weight': 'bold', 'size': lsize + 2.5, }
    tem_x = np.linspace(1993, 2020, 28)
    tem_x_filt = np.linspace(1993, 2020, 200)


    fig = plt.figure(figsize=(14, 14), dpi=200)  # 设置大小和分辨率
    G = gridspec.GridSpec(15, 2, figure=fig)
    ax = plt.subplot(G[0:4, :])

    m_medium = 161
    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_8/dpco2_current.mat')
    bin_edge = np.squeeze(tem_data['bin_edge'])
    m_hist = np.squeeze(tem_data['m_hist'])

    plt.plot(bin_edge[:m_medium], m_hist[:m_medium],  c=c2[0], linewidth=2.0, alpha=al)
    plt.plot(bin_edge[m_medium - 1:], m_hist[m_medium - 1:],  c=c2[1], linewidth=2.0, alpha=al)

    plt.fill_between(bin_edge[:m_medium], np.zeros(bin_edge[:m_medium].__len__()), m_hist[:m_medium],  color=c2[0], alpha=0.2)
    plt.fill_between(bin_edge[m_medium - 1:], np.zeros(bin_edge[m_medium - 1:].__len__()), m_hist[m_medium - 1:], color=c2[1],
                     alpha=0.2)

    plt.ylim([0, 0.035])
    plt.vlines(0, 0, 1, colors='k', linestyles='--', linewidth=1.)
    plt.xlabel(r'$DpCO_2\ (\mu atm)$',fontsize=lsize)
    plt.ylabel('Probability', fontsize=lsize)
    plt.text(0.25, 0.6, r'$Undersaturated$', fontsize=lsize, c=c2[0], transform=ax.transAxes, horizontalalignment='center')
    plt.text(0.8, 0.6, r'$Supersaturated$', fontsize=lsize, c=c2[1], transform=ax.transAxes, horizontalalignment='center')
    plt.text(0.25, 0.5, r'$20.5\%$', fontsize=lsize, c=c2[0], transform=ax.transAxes, horizontalalignment='center')
    plt.text(0.8, 0.5, r'$79.5\%$', fontsize=lsize, c=c2[1], transform=ax.transAxes, horizontalalignment='center')
    ax.tick_params(labelsize=lsize)
    plt.title('a Prestorm', x=0.01, y=0.95, font=font_l, horizontalalignment='left', verticalalignment='top')


    ax = plt.subplot(G[4:8, :])
    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_8/dpco2_cmip.mat')
    bin_edge = np.squeeze(tem_data['bin_edge'])
    dpco2_all_line_mean = np.squeeze(tem_data['dpco2_all_line_mean'])

    plt.fill_between(bin_edge, dpco2_all_line_mean[0, :], dpco2_all_line_mean[1, :], where= dpco2_all_line_mean[0, :] < dpco2_all_line_mean[1, :],
                     hatch='..', facecolor='white', alpha=0.2)
    plt.fill_between(bin_edge, dpco2_all_line_mean[0, :], dpco2_all_line_mean[1, :], where= dpco2_all_line_mean[0, :] >= dpco2_all_line_mean[1, :],
                     hatch='//', facecolor='white', alpha=0.2)

    m_hist = dpco2_all_line_mean[0, :]
    plt.plot(bin_edge[:m_medium], m_hist[:m_medium], c=c2[0], linewidth=2.0, alpha=al)
    plt.plot(bin_edge[m_medium - 1:], m_hist[m_medium - 1:], c=c2[1], linewidth=2.0, alpha=al)
    plt.fill_between(bin_edge[:m_medium], np.zeros(bin_edge[:m_medium].__len__()), m_hist[:m_medium],
                     color=c2[0], alpha=0.2)
    plt.fill_between(bin_edge[m_medium - 1:], np.zeros(bin_edge[m_medium - 1:].__len__()), m_hist[m_medium - 1:],
                     color=c2[1], alpha=0.2)

    m_hist = dpco2_all_line_mean[1, :]
    l2, = plt.plot(bin_edge[:m_medium], m_hist[:m_medium], c=c2[0], linewidth=2.0, alpha=al, linestyle='--')
    plt.plot(bin_edge[m_medium - 1:], m_hist[m_medium - 1:], c=c2[1], linewidth=2.0, alpha=al, linestyle='--')
    plt.fill_between(bin_edge[:m_medium], dpco2_all_line_mean[0, :][:m_medium], dpco2_all_line_mean[1, :][:m_medium],
                     where=dpco2_all_line_mean[0, :][:m_medium] < dpco2_all_line_mean[1, :][:m_medium],
                     color=c2[0], alpha=0.2)
    plt.fill_between(bin_edge[m_medium - 1:], dpco2_all_line_mean[0, :][m_medium - 1:], dpco2_all_line_mean[1, :][m_medium - 1:],
                     where=dpco2_all_line_mean[0, :][m_medium - 1:] < dpco2_all_line_mean[1, :][m_medium - 1:],
                     color=c2[1], alpha=0.2)

    l1, = plt.plot([0, 0], [0, 0], c='k', linestyle='-', linewidth=2.0, alpha=al)
    l2, = plt.plot([0, 0], [0, 0], c='k', linestyle='--', linewidth=2.0, alpha=al)
    l3 = plt.fill_between([0, 0], [0, 0], [0, 0], hatch='..', facecolor='#F1F1F1', alpha=0.1)
    l4 = plt.fill_between([0, 0], [0, 0], [0, 0], hatch='//', facecolor='#F1F1F1', alpha=0.1)
    plt.legend([l1, l2, l3, l4,], ['1993-2020', '2070-2100', 'Increase', 'Decrease'], handlelength=1.6, fontsize=lsize, frameon=False, loc=1)
    plt.ylim([0, 0.03])
    plt.vlines(0, 0, 1, colors='k', linestyles='--', linewidth=1.)
    plt.xlabel(r'$DpCO_2\ (\mu atm)$', fontsize=lsize)
    plt.ylabel(r'Probability', fontsize=lsize)
    ax.tick_params(labelsize=lsize)
    plt.title('b climatology', x=0.01, y=0.95, font=font_l, horizontalalignment='left', verticalalignment='top')
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.01))



    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    tem_data = scio.loadmat('./intermediate_file/Extended_Fig_8/dpco2_line.mat')
    dpco2_line_year = tem_data['dpco2_line_year']

    ax = plt.subplot(G[8:12, :])
    tem_data = np.squeeze(dpco2_line_year[:, 2])
    l2, = plt.plot(tem_x, tem_data, linewidth = 2., c='k')
    plt.plot(tem_x, tem_data, '', linewidth=2., marker='o', markersize=10, c='k')
    lower, upper, m_line, popt, cov = m_fit1d(tem_x, tem_data, tem_x_filt, e=None,  error_on=False)
    cov = np.sqrt(np.diag(cov))
    print(popt, cov)
    plt.fill_between(tem_x_filt, lower, upper, color='gray', alpha = m_al)
    plt.plot(tem_x_filt, m_line, linewidth=3.2, linestyle='--', c='tab:red', alpha=0.8)
    plt.plot(tem_x_filt, upper, linewidth=1.5, linestyle='--', c='gray', alpha=m_al)
    plt.plot(tem_x_filt, lower, linewidth=1.5, linestyle='--', c='gray', alpha=m_al)
    ax.tick_params(labelsize=lsize)
    plt.xlim([1989, 2022])
    plt.ylim([0, 30])
    plt.ylabel(r'D$pCO_2\ (\mu atm)$', fontsize=lsize)
    plt.xlabel(r'Year', fontsize=lsize)
    plt.title('c Prestorm', x=0.01, y=0.95, font=font_l, horizontalalignment='left',
              verticalalignment='top')
    plt.text(0.05, 0.08, r'$s = -0.04 \pm 0.05$', fontsize=lsize, c='k', transform=ax.transAxes)


    ax = plt.subplot(G[12:, 0])
    tem_data = dpco2_line_year[:, 0]
    l1, = plt.plot(tem_x, tem_data, linewidth = 2., c=c2[0])
    plt.plot(tem_x, tem_data, '', linewidth=2., marker='o', markersize=10, c=c2[0])
    lower, upper, m_line, popt, cov = m_fit1d(tem_x, tem_data, tem_x_filt, e=None,  error_on=False)
    cov = np.sqrt(np.diag(cov))
    print(popt, cov)
    plt.fill_between(tem_x_filt, lower, upper, color=c2[0], alpha = m_al)
    plt.plot(tem_x_filt, m_line, linewidth=3.2, linestyle='--', c=c2[0], alpha=0.8)
    plt.plot(tem_x_filt, upper, linewidth=1.5, linestyle='--', c=c2[0], alpha=m_al)
    plt.plot(tem_x_filt, lower, linewidth=1.5, linestyle='--', c=c2[0], alpha=m_al)
    ax.tick_params(labelsize=lsize)
    plt.xlim([1989, 2022])
    plt.ylim([-30, 0])
    plt.ylabel(r'D$pCO_2\ (\mu atm)$', fontsize=lsize)
    plt.xlabel(r'Year', fontsize=lsize)
    plt.title('d Prestorm', x=0.02, y=0.945, font=font_l, horizontalalignment='left',
              verticalalignment='top')
    plt.text(0.6, 0.1, r'$s = -0.09 \pm 0.05$', fontsize=lsize, c=c2[0], transform=ax.transAxes)
    plt.legend([l1], [r'$Undersaturated$', ],
               fontsize=lsize, loc=1, framealpha=0.6, edgecolor='gray'
               , frameon=False, ncol=1)


    ax = plt.subplot(G[12:, 1])
    tem_data = dpco2_line_year[:, 1]
    l2, = plt.plot(tem_x, tem_data, linewidth = 2., c=c2[1])
    plt.plot(tem_x, tem_data, '', linewidth=2., marker='o', markersize=10, c=c2[1])
    lower, upper, m_line, popt, cov = m_fit1d(tem_x, tem_data, tem_x_filt, e=None,  error_on=False)
    cov = np.sqrt(np.diag(cov))
    print(popt, cov)
    plt.fill_between(tem_x_filt, lower, upper, color=c2[1], alpha = m_al)
    plt.plot(tem_x_filt, m_line, linewidth=3.2, linestyle='--', c=c2[1], alpha=0.8)
    plt.plot(tem_x_filt, upper, linewidth=1.5, linestyle='--', c=c2[1], alpha=m_al)
    plt.plot(tem_x_filt, lower, linewidth=1.5, linestyle='--', c=c2[1], alpha=m_al)
    ax.tick_params(labelsize=lsize)
    plt.xlim([1989, 2022])
    plt.ylim([5, 35])
    plt.yticks(np.linspace(5,35,4), [5, 15, 25, 35])
    plt.xlabel(r'Year', fontsize=lsize)
    plt.title('e Prestorm', x=0.02, y=0.945, font=font_l, horizontalalignment='left',
              verticalalignment='top')
    plt.text(0.6, 0.1, r'$s = 0.03 \pm 0.04$', fontsize=lsize, c=c2[1], transform=ax.transAxes)
    plt.legend([l2], [r'$Supersaturated$', ],
               fontsize=lsize, loc=1, framealpha=0.6, edgecolor='gray',
               frameon=False, ncol=1)


    plt.tight_layout()
    plt.show()

    time_end = Time.time()
    print('cost %f s' % (time_end-time_begin))