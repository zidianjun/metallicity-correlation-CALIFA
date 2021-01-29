
from post_process import DataHub
from diagnostics import split_uf_array, merge_uf_array
import config

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import pearsonr
from npdfplot import npdfplot

def bootstrap(x_axis, times=50):
    xdata, ydata = [], []
    pr = []
    dh = DataHub()
    for i in range(times):
        x = np.random.normal(dh.col(x_axis), dh.col('error_'+x_axis))
        xdata = np.append(xdata, x)
        l = dh.rand_corr_len()
        ydata = np.append(ydata, l)
        if x_axis == 'SFR':
            mask = (x > 0) # remove nan in SFR
            pr.append(pearsonr(np.log10(x[mask]), np.log10(l[mask]))[0])
        else:
            pr.append(pearsonr(x, np.log10(l))[0])
        print i
    print np.mean(pr), np.std(pr)
    return xdata, ydata

def hist2d(x_axis, y_range=(-1.5, .9, 20)):
    xdata, ydata = bootstrap(x_axis)

    fig, ax = plt.subplots(figsize=(13, 9))
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.ylim(.04, 10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    
    if x_axis == 'SFR':
        x_range = (-1.2, 1.5, 20)
        plt.xlabel("SFR (M$_{\odot}$/yr)", fontsize=20)
    elif x_axis == 'Mass':
        x_range = (9., 12., 20)
        plt.xlabel("Stellar mass (M$_{\odot}$)", fontsize=20)
    elif x_axis == 'age':
        x_range = (7.5, 10., 20)
        plt.xlabel("Luminosity-weighted age (yr)", fontsize=20)
    else:
        x_range = (2., 4., 20)
        plt.xlabel("Surface density (M$_{\odot}$/pc$^2$)", fontsize=20)
    plt.ylabel("Correlation length (kpc)", fontsize=20)

    x0, x1, nx = x_range
    y0, y1, ny = y_range

    dx, dy = (x1 - x0) / nx, (y1 - y0) / ny

    x, y = np.meshgrid(np.arange(x0, x1, dx), np.arange(y0, y1, dy))
    x, y = np.reshape(x, [1, -1])[0], np.reshape(y, [1, -1])[0]
    count = np.zeros(nx * ny)
    for i in range(len(xdata)):
        if xdata[i] > 0. and ydata[i] > 0.:
            if x_axis == 'SFR':
                ind_x = int((np.log10(xdata[i]) - x0) / dx)
            else:
                ind_x = int((xdata[i] - x0) / dx)
            ind_y = int((np.log10(ydata[i]) - y0) / dy)
            ind = ind_y * nx + ind_x
            if 0 < ind < nx * ny:
                count[ind] += 1

    c_axis = np.log10(count/max(count))
    mask = c_axis >= -1.5
    fig = plt.scatter(10**x[mask], 10**y[mask], c=c_axis[mask],
                      s=8.2e2, edgecolors='face', marker='s', cmap=plt.cm.afmhot_r)
    cbar = plt.colorbar(fig)
    cbar.set_label('log probability density', size=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xlim(10**x0, 10**x1)

gal_prop = 'Mass'
hist2d(gal_prop)
plt.savefig(config.savefig_path + gal_prop + '.pdf')
gal_prop = 'SFR'
hist2d(gal_prop)
plt.savefig(config.savefig_path + gal_prop + '.pdf')

#plt.show()


