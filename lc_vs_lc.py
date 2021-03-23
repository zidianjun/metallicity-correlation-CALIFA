
import config

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from scipy.stats import pearsonr
from utils import open_full_chain
from post_process import DataHub

def bootstrap(suffix, dh=DataHub(), times=50):
    pr = []
    for i in range(times):
        l1, l2 = dh.rand_corr_len(), dh.rand_corr_len(suffix=suffix)
        pr.append(pearsonr(np.log10(l1), np.log10(l2))[0])
        print i
    print np.mean(pr), np.std(pr)


suffix = '_Ka03'
lcorr1 = pd.read_csv(config.output_path + '/correlation_length.csv')
lcorr2 = pd.read_csv(config.output_path + '/correlation_length' + suffix + '.csv')

dh = DataHub()

l1_50, l1_16, l1_84 = lcorr1.l_50[dh.mask], lcorr1.l_16[dh.mask], lcorr1.l_84[dh.mask]
l2_50, l2_16, l2_84 = lcorr2.l_50[dh.mask], lcorr2.l_16[dh.mask], lcorr2.l_84[dh.mask]

x = np.arange(.01, 10, .01)

plt.subplots(figsize=(10, 10))

ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
ax.vlines(l1_50, l2_16, l2_84, color='gray')
ax.hlines(l2_50, l1_16, l1_84, color='gray')
ax.scatter(l1_50, l2_50, color='k')
if suffix == '_adp':
    ax.set_xlabel('$l_{\mathrm{corr}}$ [fixed bin width] (kpc)', fontsize=20)
    ax.set_ylabel('$l_{\mathrm{corr}}$ [adaptive bin width] (kpc)', fontsize=20)
else:
    ax.set_xlabel('$l_{\mathrm{corr}}$ [Kewley line] (kpc)', fontsize=20)
    ax.set_ylabel('$l_{\mathrm{corr}}$ [Kauffmann line] (kpc)', fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.set_xlim(.06, 10)
ax.set_ylim(.06, 10)
ax.plot(x, x, color='gray', linestyle='--')


#bootstrap(suffix)
plt.savefig(config.savefig_path + suffix + '.pdf')
#plt.show()
