
import config
import constant
from post_process import DataHub
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import pearsonr

def bootstrap(d1, d2, times=20):
    xdata, ydata = [], []
    pr = []
    dh = DataHub()
    l50 = dh.corr_len()[0]
    for i in range(times):
        l1 = dh.rand_corr_len(diag=d1)
        xdata = np.append(xdata, l1)
        l2 = dh.rand_corr_len(diag=d2)
        ydata = np.append(ydata, l2)
        pr.append(pearsonr(np.log10(l1), np.log10(l2))[0])
        #print(i)
    print(np.mean(pr), np.std(pr))
    return xdata, ydata

def hist2d(xdata, ydata):
    x, y = np.meshgrid(np.arange(-1.3, 1.2, .25), np.arange(-1.3, 1.2, .25))
    x, y = np.reshape(x, [1, -1])[0], np.reshape(y, [1, -1])[0]
    count = np.zeros(100)
    for i in range(len(xdata)):
        if xdata[i] > 0 and ydata[i] > 0:
            ind_x = int((np.log10(xdata[i]) + 1.425) / .25)
            ind_y = int((np.log10(ydata[i]) + 1.425) / .25)
            ind = ind_y * 10 + ind_x
            if ind < 100:
                count[ind] += 1
    
    c_axis = np.log10(count/max(count))
    mask = c_axis >= -1.5
    fig = plt.scatter(10**x[mask], 10**y[mask], c=c_axis[mask],
                s=5e2, edgecolors='face', marker='s', cmap=plt.cm.afmhot_r)
    return fig

def corner(d1, d2, num, dh=DataHub()):
    ax = plt.subplot(220+num)
    ax.set(aspect=1.0/ax.get_data_ratio())

    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.xlim(3e-2, 20)
    plt.ylim(3e-2, 20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    if num % 2 == 1:
        plt.ylabel('$l_{\mathrm{corr}}$ [%s] (kpc)' %(d2), fontsize=20)
    if num > 2:
        plt.xlabel('$l_{\mathrm{corr}}$ [%s] (kpc)' %(d1), fontsize=20)
    
    X, Y = bootstrap(d1, d2)
    fig = hist2d(X, Y)
    
    if num == 1:
        cbar = plt.colorbar(fig)
        cbar.set_label('log probability density', size=20)
        cbar.ax.tick_params(labelsize=20)
    
    x = np.arange(-.2, 1e2)
    plt.plot(x, x, color='gray', linestyle='--')
    
    print('     ')

plt.subplots(figsize=(10, 10))

corner('PPN2', 'PPO3N2', 1)
corner('PPN2', 'K19N2O2', 3)
corner('PPO3N2', 'K19N2O2', 4)

plt.savefig(config.savefig_path + 'cbar_corner.pdf')
#plt.show()

    