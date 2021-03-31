
import config
import constant
from post_process import DataHub

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import pearsonr

def bootstrap(dh=DataHub(), times=50):
    pr = []
    x = (dh.col('PSF') / np.sqrt(np.log(256)) /
         np.sqrt((dh.col('b2a') ** 2 - config.q0 ** 2) / (1 - config.q0 ** 2)))
    for i in range(times):
        y = dh.rand_corr_len() / dh.col('dist') / constant.kpc_per_Mpc
        pr.append(pearsonr(np.log10(x), np.log10(y))[0])
        #print(i)
    print(np.mean(pr), np.std(pr))


dh = DataHub()

X0 = (dh.col('PSF') / np.sqrt(np.log(256)) / 
     np.sqrt((dh.col('b2a') ** 2 - config.q0 ** 2) / (1 - config.q0 ** 2)))
Y, Y_lower, Y_upper = dh.corr_len()

Y0 = Y / dh.col('dist') / constant.kpc_per_Mpc
Y_upper0 = Y_upper / dh.col('dist') / constant.kpc_per_Mpc
Y_lower0 = Y_lower / dh.col('dist') / constant.kpc_per_Mpc

plt.subplots(figsize=(10, 10))
ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
x = np.arange(.8, 5., .01)
ax.plot(x, x, color='gray', linestyle='--')
ax.set_xlim(.8, 5)
ax.set_ylim(.1, 30)
ax.set_xlabel('$\lambda_{\mathrm{beam}}$ (arcsec)', fontsize=20)
ax.set_ylabel('$\lambda_{\mathrm{corr}}$ (arcsec)', fontsize=20)
x_ticks = [0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0]
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_ticks)
y_ticks = [0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_ticks)
ax.tick_params(axis='both', labelsize=20)

ax.scatter(X0, Y0, color='w', edgecolor='k', s=150)
ax.vlines(X0, Y_lower0, Y_upper0, color='gray', linewidth=2)


bootstrap()

plt.savefig(config.savefig_path + 'beam_size.pdf')
#plt.show()

