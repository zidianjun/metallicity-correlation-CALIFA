
import config
from post_process import DataHub

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

dh = DataHub()
HW50M, HW30M = dh.corr_scale()
corr_len = dh.corr_len()[0]

plt.subplots(figsize=(10, 10))
ax = plt.subplot(111)
ax.set_xscale("log")

plt.hist(HW30M, bins=10 ** np.arange(2, 4, .2),
         histtype='step', facecolor='m', fill=True,
         linewidth=2, edgecolor='k', alpha=.75, label='30% correlation scale\n(this work)')
plt.hist(HW50M, bins=10 ** np.arange(2, 4, .2),
         histtype='step', facecolor='c', fill=True,
         linewidth=2, edgecolor='k', alpha=.75, label='50% correlation scale\n(this work)')
plt.hist(corr_len * 1e3, bins=10 ** np.arange(2, 5, .2),
         histtype='step', facecolor='gray', hatch='\\', fill=True,
         linewidth=2, edgecolor='k', alpha=.5, label='correlation length\n(this work)')


HW50M_K20 = np.array([270, 380, 290, 230, 370, 290, 340, 370])
HW30M_K20 = np.array([510, 710, 510, 360, 760, 550, 740, 820])
plt.scatter(HW50M_K20, np.ones(8) * 45, marker='|', s=1000, color='c',
	        label='50% correlation scale\n(Kreckel et al. 2020)')
plt.scatter(HW50M_K20 - 2, np.ones(8) * 45, marker='|', s=1000, color='c',)
plt.scatter(HW50M_K20 + 2, np.ones(8) * 45, marker='|', s=1000, color='c',)
plt.scatter(HW30M_K20, np.ones(8) * 40, marker='|', s=1000, color='m',
	        label='30% correlation scale\n(Kreckel et al. 2020)')
plt.scatter(HW30M_K20 - 2, np.ones(8) * 40, marker='|', s=1000, color='m',)
plt.scatter(HW30M_K20 + 2, np.ones(8) * 40, marker='|', s=1000, color='m',)

x_ticks = [100, 200, 400, 600, 1000, 2000, 4000, 6000]
plt.xticks(x_ticks, x_ticks, fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Correlation scale & length (pc)', fontsize=20)
plt.ylabel('N', fontsize=20)
plt.xlim(100, 6000)
plt.ylim(.5, 50)
plt.legend(loc='upper right')


plt.savefig(config.savefig_path + '/corr_scale.pdf')
#plt.show()


