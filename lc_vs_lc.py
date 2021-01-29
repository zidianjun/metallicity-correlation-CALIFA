
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


suffix = '_adp'
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
plt.vlines(l1_50, l2_16, l2_84, color='gray')
plt.hlines(l2_50, l1_16, l1_84, color='gray')
plt.scatter(l1_50, l2_50, color='k')
plt.xlabel('$l_{\mathrm{corr}}$ [Kewley line] (kpc)', fontsize=20)
plt.ylabel('$l_{\mathrm{corr}}$ [Kauffmann line] (kpc)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(.06, 10)
plt.ylim(.06, 10)
plt.plot(x, x, color='gray', linestyle='--')
plt.savefig(config.savefig_path + suffix + '.pdf')


bootstrap(suffix)
#plt.show()
