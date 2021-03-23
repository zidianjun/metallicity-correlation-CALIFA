
from post_process import DataHub
from utils import read_CALIFA_catalog
import config

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import pearsonr

decom_file = read_CALIFA_catalog(name='decom.csv')
dh = DataHub()

b2a_MA17 = np.array(decom_file.b2a)
PA_MA17 = np.array(decom_file.PA)
b2a_DR3 = dh.col('b2a')
PA_DR3 = dh.col('PA') - 90


plt.figure(figsize=(12, 10))
plt.subplots_adjust(left=.12, bottom=.10, right=.98, top=.98, wspace=.3, hspace=.2)

ax = plt.subplot(221)
ax.scatter(b2a_MA17, b2a_DR3, color='w', edgecolor='k', s=50)
x = np.linspace(0, 1, 100)
ax.plot(x, x, color='gray', linestyle='--')
ax.plot(x, x-.054*2, color='gray', linestyle='--')
ax.plot(x, x+.054*2, color='gray', linestyle='--')
ax.set_xlabel("$b/a$ (MA17)", fontsize=20)
ax.set_ylabel("$b/a$ (DR3)", fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.set_xlim(.25, 1.)
ax.set_ylim(.25, 1.)

ax = plt.subplot(222)
d_b2a = b2a_MA17 - b2a_DR3
ax.hist(d_b2a[~np.isnan(d_b2a)], bins=np.arange(-.21, .23, .02),
        histtype='step', facecolor='w', fill=True, linewidth=1, edgecolor='k')
ax.set_xlabel("$\Delta$ $b/a$", fontsize=20)
ax.set_ylabel("N", fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.set_xlim(-.25, .25)

ax = plt.subplot(223)
ax.scatter(PA_MA17, PA_DR3, color='w', edgecolor='k', s=50)
x = np.linspace(-20, 200, 100)
ax.plot(x, x, color='gray', linestyle='--')
ax.plot(x, x-6.7*2, color='gray', linestyle='--')
ax.plot(x, x+6.7*2, color='gray', linestyle='--')
ax.set_xlabel("PA (MA17)", fontsize=20)
ax.set_ylabel("PA (DR3)", fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.set_xlim(-20, 200)
ax.set_ylim(-20, 200)

ax = plt.subplot(224)
d_PA = PA_MA17 - PA_DR3
ax.hist(d_PA[~np.isnan(d_PA)], bins=np.arange(-9, 11, 2),
        histtype='step', facecolor='w', fill=True, linewidth=1, edgecolor='k')
ax.set_xlabel("$\Delta$ PA", fontsize=20)
ax.set_ylabel("N", fontsize=20)
ax.tick_params(axis='both', labelsize=20)

plt.savefig(config.savefig_path + 'decom.pdf')
#plt.show()

