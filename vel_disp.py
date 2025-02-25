
from utils import read_CALIFA_catalog
import config

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

CALIFA_data = read_CALIFA_catalog()

corr_len = read_CALIFA_catalog(name='correlation_length.csv',
           path=config.output_path)

CO_data = read_CALIFA_catalog(
          name='EDGE_CO_VelDisp_Gaussian_BeamSmearingCorrected_Masked.csv')

name_list = list(corr_len.name)

x_50, x_16, x_84 = [], [], []
y_50, y_16, y_84 = [], [], []
c = []

for name in CO_data.Name:
    if (name in name_list and
        CALIFA_data[CALIFA_data.name == name].values[0, 3] > .316 and
        corr_len[corr_len.name == name].values[0, 1] > 0):

        vel_disp, e_vel_disp = CO_data[CO_data.Name == name].values[0, 1:3]
        x_50.append(vel_disp)
        x_16.append(vel_disp - e_vel_disp)
        x_84.append(vel_disp + e_vel_disp)

        y_50.append(corr_len[corr_len.name == name].values[0, 1])
        y_16.append(corr_len[corr_len.name == name].values[0, 2])
        y_84.append(corr_len[corr_len.name == name].values[0, 3])


plt.subplots(figsize=(10, 10))
ax = plt.subplot(111)
ax.set_yscale("log")
ax.scatter(x_50, y_50, color='w', edgecolor='k', s=100)
ax.vlines(x_50, y_16, y_84, color='gray')
ax.hlines(y_50, x_16, x_84, color='gray')

x = np.arange(0, 40, .1)
ax.plot(x, np.sqrt(1e-3/3*400*x*3), color='k', linestyle='--',
         label='scale height = 400pc, SF duration = 3Gyr')
ax.plot(x, np.sqrt(1e-3/3*180*x*2), color='k', linestyle='-.',
         label='scale height = 180pc, SF duration = 2Gyr')
ax.plot(x, np.sqrt(1e-3/3*80*x*1), color='k', linestyle=':',
         label='scale height = 80pc,   SF duration = 1Gyr')
ax.legend(loc='upper left', prop={'size': 15})


ax.set_xlim(0, 40)
ax.set_ylim(0.1, 10)

ax.set_xlabel('CO velocity dispersion (km s$^{-1}$)', fontsize=20)
ax.set_ylabel('Correlation length (kpc)', fontsize=20)
ax.tick_params(axis='both', labelsize=20)


plt.savefig(config.savefig_path + 'vel_disp.pdf')
#plt.show()







