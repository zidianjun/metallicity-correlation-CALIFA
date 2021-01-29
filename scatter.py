
from post_process import DataHub
import config
from utils import read_CALIFA_catalog

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def scatter_plot(dh=DataHub()):
    plt.subplots(figsize=(10, 8))
    ax = plt.subplot(111)
    plt.xlim(0, 160)
    plt.ylim(0, 45)
    g = read_CALIFA_catalog()
    x0, y0 = g.dist, g.R
    fig = plt.scatter(x0, y0, s=75, color='w', edgecolor='k')
    x1 = dh.col('dist')
    y1 = dh.col('R')
    fig = plt.scatter(x1, y1, s=75, color='mediumblue', edgecolor='k')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('distance (Mpc)', fontsize=20)
    plt.ylabel('$\mathrm{R}_{\mathrm{iso, r}}$ (kpc)', fontsize=20)
    plt.savefig(config.savefig_path + 'Re.pdf')

def bootstrap(dh=DataHub(), times=50):
    mat = []
    for i in range(times):
        mat.append(dh.rand_corr_len())
        print i
    return np.array(mat)

def mean_std(mat, mask):
    m = []
    for i in range(len(mat)):
        m.append(np.mean(mat[i][mask]))
    return np.mean(m), np.std(m)

def T_type(mat, dh=DataHub()):
    plt.subplots(figsize=(10, 5))
    plt.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85)

    n_type = 8
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('T-type') == (10 + i))
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)

    plt.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    plt.vlines(x, y - e, y + e)
    plt.xlim(0, 9)
    plt.ylim(.0, 2.0)
    plt.xlabel('Hubble type', fontsize=20)
    plt.ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    labels = ['', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm', '']
    plt.xticks(ticks, labels, fontsize=20)
    plt.yticks(fontsize=20)
    for i in range(n_type):
        plt.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)
    plt.savefig(config.savefig_path + '/morph/hubble_type.pdf')


def bar(mat, dh=DataHub()):
    plt.subplots(figsize=(10, 5))
    plt.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85)
    
    n_type = 3
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('bar') == i)
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)
    
    plt.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    plt.vlines(x, y - e, y + e)
    plt.xlim(0, 3)
    plt.ylim(.5, 1.5)
    plt.xlabel('barredness', fontsize=20)
    plt.ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    labels = ['', '', 'A', '', 'AB', '', 'B', '', '']
    plt.xticks(ticks, labels, fontsize=20)
    plt.yticks(fontsize=20)
    for i in range(n_type):
        plt.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)
    plt.savefig(config.savefig_path + '/morph/bar.pdf')

def merging(mat, dh=DataHub()):
    plt.subplots(figsize=(10, 5))
    plt.subplots_adjust(left=.15, bottom=.15, right=.85, top=.85)
    
    n_type = 2
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('merge') == i)
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)
    
    plt.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    plt.vlines(x, y - e, y + e)
    plt.xlim(0, 3)
    plt.ylim(.5, 3.)
    plt.xlabel('interaction state', fontsize=20)
    plt.ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
    labels = ['', '', 'without merging', '', 'with merging', '', '']
    plt.xticks(ticks, labels, fontsize=20)
    plt.yticks(fontsize=20)
    for i in range(n_type):
        plt.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)
    plt.savefig(config.savefig_path + '/morph/merging.pdf')

#scatter_plot()


mat = bootstrap()
T_type(mat)
bar(mat)
merging(mat)

#plt.show()


