
import config
import constant
from post_process import DataHub

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import pearsonr

def bootstrap(dh=DataHub(), times=50):
    pr = []
    x = (dh.col('FWHM') / np.sqrt(np.log(256)) /
         np.sqrt((dh.col('b2a') ** 2 - config.q0 ** 2) / (1 - config.q0 ** 2)))
    for i in range(times):
        y = dh.rand_corr_len() / dh.col('dist') / constant.kpc_per_Mpc
        pr.append(pearsonr(np.log10(x), np.log10(y))[0])
        print i
    print np.mean(pr)
    print np.std(pr)


dh = DataHub()

X0 = (dh.col('FWHM') / np.sqrt(np.log(256)) / 
     np.sqrt((dh.col('b2a') ** 2 - config.q0 ** 2) / (1 - config.q0 ** 2)))
Y, Y_lower, Y_upper = dh.corr_len()

Y0 = Y / dh.col('dist') / constant.kpc_per_Mpc
Y_upper0 = Y_upper / dh.col('dist') / constant.kpc_per_Mpc
Y_lower0 = Y_lower / dh.col('dist') / constant.kpc_per_Mpc

plt.subplots(figsize=(10, 10))
ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
x = np.arange(.1, 10)
plt.plot(x, x, color='gray', linestyle='--')
plt.xlim(.8, 5)
plt.ylim(.1, 30)
plt.xlabel('$\lambda_{\mathrm{beam}}$ (arcsec)', fontsize=20)
plt.ylabel('$\lambda_{\mathrm{corr}}$ (arcsec)', fontsize=20)
x_ticks = [0.8, 0.9, 1.0, 2.0, 3.0, 4.0]
plt.xticks(x_ticks, x_ticks, fontsize=20)
y_ticks = [0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0]
plt.yticks(y_ticks, y_ticks, fontsize=20)
plt.yticks(fontsize=20)

plt.vlines(X0, Y_lower0, Y_upper0, color='k', linewidth=2)
plt.scatter(X0, Y0, color='w', edgecolor='k', s=150)
print pearsonr(np.log10(X0), np.log10(Y0))[0]

bootstrap()

plt.savefig(config.savefig_path + 'beam_size.pdf')
#plt.show()

