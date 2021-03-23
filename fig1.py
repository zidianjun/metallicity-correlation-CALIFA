
from post_process import DataHub
import config
from utils import read_CALIFA_catalog

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def panel1(dh=DataHub()):
    plt.figure(figsize=(10, 10))
    ax = plt.axes([.15, .15, .8, .8])
    plt.xlim(0, 160)
    plt.ylim(0, 23)
    g = read_CALIFA_catalog()
    x0, y0 = g.dist, g.R_25
    fig = plt.scatter(x0, y0, s=75, color='w', edgecolor='k')
    x1 = dh.col('dist')
    y1 = dh.col('R_25')
    fig = plt.scatter(x1, y1, s=75, color='royalblue', edgecolor='k')
    plt.scatter([9.77, 14.4, 11.9, 10.1, 10.6, 16.8, 15.8, 9.95],
                [14.2, 6.3, 10.7, 9.4, 15.7, 12.2, 18.8, 8.7],
                s=150, marker='D', color='tomato', edgecolor='k')
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('distance (Mpc)', fontsize=30)
    plt.ylabel('$\mathrm{R}_{25}$ (kpc)', fontsize=30)
    
    plt.savefig(config.savefig_path + 'R_25.pdf')

def panel2(dh=DataHub()):
    g = read_CALIFA_catalog()
    x0, y0 = 10 ** g.Mass, 10 ** g.lSFR

    x, y = 10 ** dh.col('Mass'), 10 ** dh.col('lSFR')

    left, width = 0.2, 0.55
    bottom, height = 0.15, 0.6
    spacing = 0.

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    plt.figure(figsize=(11, 10))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.set_xscale("log")
    ax_scatter.set_yscale("log")
    ax_scatter.tick_params(top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(labelleft=False)

    ax_scatter.scatter(x0, y0, s=50, color='w', edgecolor='k')
    ax_scatter.scatter(x, y, s=50, color='royalblue', edgecolor='k')
    ax_scatter.set_xlim((2e8, 2e12))
    ax_scatter.set_ylim((1e-3, 4e1))
    ax_scatter.set_xlabel("Stellar mass (M$_{\odot}$)", fontsize=30)
    ax_scatter.set_ylabel("SFR (M$_{\odot}$ yr$^{-1}$)", fontsize=30)
    ax_scatter.tick_params(axis='both', labelsize=30)

    ax_histx.hist(x.tolist() * 153, bins=10 ** np.arange(8.3, 12., .34),
                  histtype='step', fill=True,
                  facecolor='royalblue', edgecolor='k', linewidth=1.)
    ax_histx.hist(x0.tolist() * 33, bins=10 ** np.arange(8.3, 12., .34),
                  histtype='step', fill=False,
                  edgecolor='k', linewidth=1.5)
    ax_histx.set_xscale("log")
    ax_histx.set_ylim(0, 6000)
    ticks = [.5*33*153, 33*153]
    ax_histx.set_yticks(ticks)
    ax_histx.set_yticklabels([0.5, 1.0])
    ax_histx.set_ylabel('Normalised\ndistribution', fontsize=30)
    ax_histx.tick_params(axis='y', labelsize=30)

    ax_histy.hist(y.tolist() * 142, bins=10 ** np.arange(-3, 2.5, .5),
                  orientation='horizontal', histtype='step', fill=True,
                  facecolor='royalblue', edgecolor='k', linewidth=1.)
    ax_histy.hist(y0.tolist() * 39, bins=10 ** np.arange(-3, 2.5, .5),
                  orientation='horizontal', histtype='step', fill=False,
                  edgecolor='k', linewidth=1.5)
    ax_histy.set_yscale("log")
    ax_histy.set_xlim(0, 6000)
    ticks = [.5*39*142, 39*142]
    ax_histy.set_xticks(ticks)
    ax_histy.set_xticklabels([0.5, 1.0])
    ax_histy.set_xlabel('Normalised\ndistribution', fontsize=30)
    ax_histy.tick_params(axis='x', labelsize=30)

    plt.savefig(config.savefig_path + '/hist_Mass_SFR.pdf')



def panel3(dh=DataHub()):
    g = read_CALIFA_catalog()
    x0, y0 = 10 ** g.Mass, g.R_25


    x, y = 10 ** dh.col('Mass'), dh.col('R_25')

    left, width = 0.2, 0.55
    bottom, height = 0.15, 0.6
    spacing = 0.

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    plt.figure(figsize=(10, 10))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.set_xscale("log")
    ax_scatter.set_yscale("log")
    ax_scatter.tick_params(top=True, right=True)
    y_ticks = [2, 4, 6, 8, 10, 20]
    ax_scatter.set_yticks(y_ticks)
    ax_scatter.set_yticklabels(y_ticks)

    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(labelleft=False)

    ax_scatter.scatter(x0, y0, s=50, color='w', edgecolor='k')
    ax_scatter.scatter(x, y, s=50, color='royalblue', edgecolor='k')
    ax_scatter.scatter(10 ** np.array([10.2, 9.8, 10.2, 9.6, 10.6, 10.5, 10.4, 9.6]),
                                      [14.2, 6.3, 10.7, 9.4, 15.7, 12.2, 18.8, 8.7],
                       s=100, marker='D', color='tomato', edgecolor='k')  # Kreckel et al. (2020)
    ax_scatter.set_xlim((2e8, 7e11))
    ax_scatter.set_ylim((2, 25))
    ax_scatter.set_xlabel("Stellar mass (M$_{\odot}$)", fontsize=30)
    ax_scatter.set_ylabel("$\mathrm{R}_{25}$ (kpc)", fontsize=30)
    ax_scatter.tick_params(axis='both', labelsize=30)

    ax_histx.hist(x.tolist() * 153, bins=10 ** np.arange(8.3, 12., .34),
                  histtype='step', fill=True,
                  facecolor='royalblue', edgecolor='k', linewidth=1.)
    ax_histx.hist(x0.tolist() * 33, bins=10 ** np.arange(8.3, 12., .34),
                  histtype='step', fill=False,
                  edgecolor='k', linewidth=1.5)
    ax_histx.set_xscale("log")
    ax_histx.set_ylim(0, 6000)
    ticks = [.5*33*153, 33*153]
    ax_histx.set_yticks(ticks)
    ax_histx.set_yticklabels([0.5, 1.0])
    ax_histx.set_ylabel('Normalised\ndistribution', fontsize=30)
    ax_histx.tick_params(axis='y', labelsize=30)

    ax_histy.hist(y.tolist() * 168, bins=10 ** np.arange(.3, 1.57, .17),
                  orientation='horizontal', histtype='step', fill=True,
                  facecolor='royalblue', edgecolor='k', linewidth=1.)
    ax_histy.hist(y0.tolist() * 40, bins=10 ** np.arange(.3, 1.57, .17),
                  orientation='horizontal', histtype='step', fill=False,
                  edgecolor='k', linewidth=1.5)
    ax_histy.set_yscale("log")
    ax_histy.set_xlim(0, 7000)
    ticks = [.5*40*168, 40*168]
    ax_histy.set_xticks(ticks)
    ax_histy.set_xticklabels([0.5, 1.0])
    ax_histy.set_xlabel('Normalised\ndistribution', fontsize=30)
    ax_histy.tick_params(axis='x', labelsize=30)

    plt.savefig(config.savefig_path + '/hist_Mass_R.pdf')



panel1()
panel2()
panel3()
#plt.show()




