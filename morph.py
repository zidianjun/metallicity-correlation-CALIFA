
from post_process import DataHub
import config
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')


def bootstrap(dh=DataHub(), times=50):
    mat = []
    for i in range(times):
        mat.append(dh.rand_corr_len())
        print(i)
    return np.array(mat)

def mean_std(mat, mask):
    m = []
    for i in range(len(mat)):
        m.append(np.mean(mat[i][mask]))
    return np.mean(m), np.std(m)

def T_type(ax, mat, dh=DataHub()):

    n_type = 8
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('T-type') == (10 + i))
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)

    ax.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    ax.vlines(x, y - e, y + e)
    ax.set_xlim(0, 9)
    ax.set_ylim(.0, 2.0)
    ax.set_xlabel('Hubble type', fontsize=20)
    ax.set_ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    labels = ['', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm', '']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='both', labelsize=20)
    for i in range(n_type):
        ax.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)


def bar(ax, mat, dh=DataHub()):
    
    n_type = 3
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('bar') == i)
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)
    
    ax.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    ax.vlines(x, y - e, y + e)
    ax.set_xlim(0, 3)
    ax.set_ylim(.5, 1.5)
    ax.set_xlabel('barredness', fontsize=20)
    ax.set_ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    labels = ['', '', 'A', '', 'AB', '', 'B', '', '']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='both', labelsize=20)
    for i in range(n_type):
        ax.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)


def merging(ax, mat, dh=DataHub()):
    
    n_type = 2
    x, y, e, n = np.zeros(n_type), np.zeros(n_type), np.zeros(n_type), np.zeros(n_type)

    for i in range(n_type):
        x[i] = 1 + i
        mask = (dh.col('merge') == i)
        y[i], e[i] = mean_std(mat, mask)
        n[i] = np.sum(mask)
    
    ax.scatter(x, y, s=50, marker='o', color='w', edgecolor='k')
    ax.vlines(x, y - e, y + e)
    ax.set_xlim(0, 3)
    ax.set_ylim(.5, 3.)
    ax.set_xlabel('interaction state', fontsize=20)
    ax.set_ylabel('correlation length (kpc)', fontsize=20)
    ticks = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
    labels = ['', '', 'without merging', '', 'with merging', '', '']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='both', labelsize=20)
    for i in range(n_type):
        ax.annotate('%d' %(n[i]), xy=(x[i] + .1, y[i]), xytext=(x[i] + .1, y[i]), fontsize=15)


scatter_plot()

'''
mat = bootstrap(times=2)

plt.figure(figsize=(10, 16))
plt.subplots_adjust(left=.15, bottom=.08, right=.95, top=.98, hspace=.3)

ax = plt.subplot(311)
T_type(ax, mat)
ax = plt.subplot(312)
bar(ax, mat)
ax = plt.subplot(313)
merging(ax, mat)
'''



'''
dh = DataHub()
bar = read_CALIFA_catalog(name='BAR.csv')
l_50, l_16, l_84 = dh.corr_len()

plt.figure(figsize=(16, 6))
plt.subplots_adjust(left=.05, bottom=.10, right=.95, top=.98, wspace=.0)

ax = plt.subplot(131)
ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(bar.Rbarr * 4.848e-3 * dh.col('dist'), l_50)
ax.vlines(bar.Rbarr * 4.848e-3 * dh.col('dist'), l_16, l_84)
ax.hlines(l_50, (bar.Rbarr - bar.e_Rbarr) * 4.848e-3 * dh.col('dist'),
                (bar.Rbarr + bar.e_Rbarr) * 4.848e-3 * dh.col('dist'))
ax.set_ylabel('Correlation length')
ax.set_xlabel('R-band bar radius (kpc)')
ax.set_xlim(1, 40)
ax.set_ylim(.03, 10)

ax = plt.subplot(132)
ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(bar.babarr, l_50)
ax.vlines(bar.babarr, l_16, l_84)
ax.hlines(l_50, (bar.babarr - bar.e_babarr), (bar.babarr + bar.e_babarr))
ax.set_xlabel('R-band bar axis ratio')
ax.set_yticks([])
ax.set_xlim(5e-2, 2)
ax.set_ylim(.03, 10)

ax = plt.subplot(133)
ax.set_xscale('log')
ax.set_yscale('log')
ax.scatter(bar.BarTr, l_50)
ax.vlines(bar.BarTr, l_16, l_84)
ax.set_xlabel('R-band Bar/Total luminosity ratio')
ax.set_yticks([])
ax.set_xlim(-.05, .55)
ax.set_ylim(.03, 10)


plt.savefig(config.savefig_path + 'test.pdf')
#plt.show()
'''