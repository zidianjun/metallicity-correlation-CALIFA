
import config
from post_process import DataHub

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')


def corr_scale_plot(dh=DataHub()):

    plt.figure(figsize=(20, 8))
    plt.subplots_adjust(left=.05, bottom=.12, right=.93, top=.94, wspace=.2)

    HW50M, HW30M = dh.corr_scale()
    corr_len = dh.corr_len()[0]

    mask = (dh.col('Mass') > 9.6) & (dh.col('Mass') < 10.6)

    ax = plt.subplot(121)
    ax.set_xscale("log")

    ax.hist(HW30M, bins=10 ** np.arange(2, 4, .2),
             histtype='step', facecolor='m', fill=True,
             linewidth=2, edgecolor='k', alpha=.75, label='30% correlation scale\n(this work)')
    ax.hist(HW50M, bins=10 ** np.arange(2, 4, .2),
             histtype='step', facecolor='c', fill=True,
             linewidth=2, edgecolor='k', alpha=.75, label='50% correlation scale\n(this work)')
    ax.hist(corr_len * 1e3, bins=10 ** np.arange(2, 5, .2),
             histtype='step', facecolor='gray', hatch='\\', fill=True,
             linewidth=2, edgecolor='k', alpha=.5, label='correlation length\n(this work)')

    HW50M_K20 = np.array([270, 380, 290, 230, 370, 290, 340, 370])
    HW30M_K20 = np.array([510, 710, 510, 360, 760, 550, 740, 820])
    ax.scatter(HW50M_K20, np.ones(8) * 45, marker='|', s=1000, color='c',
                label='50% correlation scale\n(Kreckel et al. 2020)')
    ax.scatter(HW50M_K20 - 2, np.ones(8) * 45, marker='|', s=1000, color='c',)
    ax.scatter(HW50M_K20 + 2, np.ones(8) * 45, marker='|', s=1000, color='c',)
    ax.scatter(HW30M_K20, np.ones(8) * 40, marker='|', s=1000, color='m',
                label='30% correlation scale\n(Kreckel et al. 2020)')
    ax.scatter(HW30M_K20 - 2, np.ones(8) * 40, marker='|', s=1000, color='m',)
    ax.scatter(HW30M_K20 + 2, np.ones(8) * 40, marker='|', s=1000, color='m',)

    x_ticks = [100, 200, 400, 600, 1000, 2000, 4000, 6000]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_ticks)
    ax.tick_params(axis='both', labelsize=20)
    ax.set_xlabel('Correlation scale & length (pc)', fontsize=20)
    ax.set_ylabel('N', fontsize=20)
    ax.set_xlim(100, 6000)
    ax.set_ylim(.5, 50)
    ax.legend(loc='upper right', prop={'size': 15})


    ax = plt.subplot(122)
    ax.set_xscale("log")
    ax.set_yscale("log")
    dh = DataHub()

    ax.scatter(10 ** np.array([10.2, 9.8, 10.2, 9.6, 10.6, 10.5, 10.4, 9.6]),
                              [14.2, 6.3, 10.7, 9.4, 15.7, 12.2, 18.8, 8.7],
               c=np.log10([360, 710, 760, 510, 820, 740, 550, 510]),
               s=200, marker='D', edgecolor='k', cmap=plt.cm.get_cmap('RdYlBu_r'))
    fig = ax.scatter(10 ** dh.col('Mass'), dh.col('R_25'), c=np.log10(dh.cs.HW30M[dh.mask]),
                      s=100, edgecolor='k', cmap=plt.cm.get_cmap('RdYlBu_r'))
    cbar = plt.colorbar(fig, ticks=np.log10([200, 400, 600, 800, 1000]))
    cbar.ax.tick_params(labelsize=25)
    cbar.ax.set_yticklabels([200, 400, 600, 800, 1000])
    cbar.set_label('30% correlation scale (pc)', size=25)
    ax.set_xlabel("Stellar mass (M$_{\odot}$)", fontsize=25)
    ax.set_ylabel("$\mathrm{R}_{25}$ (kpc)", fontsize=25)
    ax.tick_params(axis='both', labelsize=25)
    y_ticks = [4, 6, 8, 10, 20]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_ticks)
    ax.set_xlim((1e9, 7e11))
    ax.set_ylim((3.2, 25))


corr_scale_plot()
plt.savefig(config.savefig_path + '/corr_scale_cmap.pdf')
#plt.show()


