
from post_process import DataHub
import config

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.switch_backend('agg')
from scipy.stats import pearsonr



def bootstrap(x_axis, diag='K19N2O2', err_col=True, r_norm=False, times=50):
    xdata, ydata = [], []
    pr = []
    dh = DataHub(diag=diag)
    for i in range(times):
        if err_col:
            if r_norm:
                x = (10 ** np.random.normal(dh.col(x_axis), dh.col('error_'+x_axis)) /
                         dh.col('R_25') ** 2 / np.pi)
            else:
                x = 10 ** np.random.normal(dh.col(x_axis), dh.col('error_'+x_axis))
        else:
            x = dh.col(x_axis) if x_axis == 'R_25' else 10 ** dh.col(x_axis)
        xdata = np.append(xdata, x)
        l = dh.rand_corr_len()
        ydata = np.append(ydata, l)
        pr.append(pearsonr(np.log10(x), np.log10(l))[0])
        #print(i)
    print(x_axis, np.mean(pr), np.std(pr))
    return xdata, ydata, np.mean(pr), np.std(pr)

def hist2d(ax, xdata, ydata, x_range, y_range=(-1.4, 1, 20)):
    ax.set_xscale("log")
    ax.set_yscale("log")

    x0, x1, nx = x_range
    y0, y1, ny = y_range

    dx, dy = (x1 - x0) / nx, (y1 - y0) / ny

    x, y = np.meshgrid(np.arange(x0, x1, dx), np.arange(y0, y1, dy))
    x, y = np.reshape(x, [1, -1])[0], np.reshape(y, [1, -1])[0]
    count = np.zeros(nx * ny)
    for i in range(len(xdata)):
        if xdata[i] > 0. and ydata[i] > 0.:
            ind_x = int((np.log10(xdata[i]) - x0) / dx)
            ind_y = int((np.log10(ydata[i]) - y0) / dy)
            ind = ind_y * nx + ind_x
            if 0 < ind < nx * ny:
                count[ind] += 1

    c_axis = np.log10(count/max(count))
    mask = c_axis >= -1.5
    fig = ax.scatter(10**x[mask], 10**y[mask], c=c_axis[mask],
                      s=2.8e2, edgecolors='face', marker='s', cmap=plt.cm.afmhot_r)
    #cbar = ax.colorbar(fig)
    #cbar.set_label('log probability density', size=25)
    #cbar.ax.tick_params(labelsize=25)



def ebar(log_e_arr_low, log_e_arr_upp, c, perc):
    return c * (10 ** np.abs([[np.percentile(log_e_arr_low, perc)],
                [np.percentile(log_e_arr_upp, perc)]]) - 1)

def panel(num, x_axis, x_range, c_coords, dh=DataHub(),
          err_col=True, r_norm=False):
    ax = plt.subplot(230+num)

    l_50, l_16, l_84 = dh.corr_len()
    err_prop_lcorr_low = np.log10(l_50 / l_16)
    err_prop_lcorr_upp = np.log10(l_84 / l_50)

    if x_axis == 'R_25':
        log_e_arr = np.zeros(len(dh.col(x_axis)))
    else:
        log_e_arr = dh.col('error_' + x_axis)

    xdata, ydata, pr, e_pr = bootstrap(x_axis, err_col=err_col, r_norm=r_norm)
    if x_axis == 'Mass' and r_norm: xdata /= 1e6
    hist2d(ax, xdata, ydata, (x_range[0], x_range[1], 20))
    cx, cy = c_coords
    ex_16 = ebar(log_e_arr, log_e_arr, cx, 16)
    ex_50 = ebar(log_e_arr, log_e_arr, cx, 50)
    ex_84 = ebar(log_e_arr, log_e_arr, cx, 84)
    ey_16 = ebar(err_prop_lcorr_low, err_prop_lcorr_upp, cy, 16)
    ey_50 = ebar(err_prop_lcorr_low, err_prop_lcorr_upp, cy, 50)
    ey_84 = ebar(err_prop_lcorr_low, err_prop_lcorr_upp, cy, 84)
    ax.errorbar(np.array([cx]), np.array([cy]), xerr=ex_16, yerr=ey_16, ecolor='gray', capsize=4)
    ax.errorbar(np.array([cx]), np.array([cy]), xerr=ex_50, yerr=ey_50, ecolor='gray', capsize=4)
    ax.errorbar(np.array([cx]), np.array([cy]), xerr=ex_84, yerr=ey_84, ecolor='gray', capsize=4)
    ax.set_ylabel("Correlation length (kpc)" * (num%3==1), fontsize=25)
    #ax.annotate('r=%.2f$\pm$%.2f' %(pr, e_pr), xy=(.05, .9), xytext=(.05, .9),
    #            xycoords='axes fraction', fontsize=25)
    ax.tick_params(axis='both', labelsize=25, labelleft=(num%3==1))
    if x_axis == 'R_25':
        x_ticks = [4, 6, 8, 10, 20]
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_ticks)
    ax.set_ylim(.03, 10)
    
    return ax

def hist2d_Mass_SFR(diag='K19N2O2', dh=DataHub()):
    plt.figure(figsize=(16, 12))
    plt.subplots_adjust(left=.10, bottom=.08, right=.98, top=.97, wspace=.0, hspace=.27)
    
    ax = panel(1, 'Mass', x_range=(9, 11.7), c_coords=(3e11, .1), dh=dh)
    ax.set_xlim(1e9, 5e11)
    ax.set_xlabel("Stellar mass (M$_{\odot}$)", fontsize=25)
    
    ax = panel(2, 'lSFR', x_range=(-1.2, 1.6), c_coords=(20, .1), dh=dh)
    ax.set_xlim(6e-2, 40)
    ax.set_xlabel("SFR (M$_{\odot}$ yr$^{-1}$)", fontsize=25)

    ax = panel(3, 'R_25', x_range=(.5, 1.4), c_coords=(20, .1), dh=dh, err_col=False)
    ax.set_xlim(3.2, 25)
    ax.set_xlabel("$\mathrm{R}_{25}$", fontsize=25)
    
    ax = panel(4, 'Mass', x_range=(1., 3.2), c_coords=(1e3, .1), dh=dh, r_norm=True)
    ax.set_xlim(12, 2e3)
    ax.set_xlabel("$\Sigma_{\mathrm{star}}$ (M$_{\odot}$ pc$^{-2}$)", fontsize=25)

    ax = panel(5, 'lSFR', x_range=(-3.2, -1), c_coords=(7e-2, .1), dh=dh, r_norm=True)
    ax.set_xlim(6e-4, 1e-1)
    ax.set_xlabel("$\Sigma_{\mathrm{SFR}}$ (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$)", fontsize=25)
       


hist2d_Mass_SFR()
plt.savefig(config.savefig_path + 'Mass_SFR.pdf')
#plt.show()
