
from utils import read_CALIFA_catalog
import config

import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

CALIFA_data = read_CALIFA_catalog()

corr_len = read_CALIFA_catalog(name='correlation_length.csv',
           path=config.output_path)

CO_data = read_CALIFA_catalog(
          name='EDGE_CO_VelDisp_Gaussian_BeamSmearingCorrected_Masked.csv',
          path=config.obj_path + '/CALIFA/')

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

        c.append((CALIFA_data[CALIFA_data.name == name].values[0, 14])) # 7: Mass, 9: SFR

plt.subplots(figsize=(10, 10))
ax = plt.subplot(111)
ax.set_yscale("log")
plt.scatter(x_50, y_50, color='w', edgecolor='k', s=100)
#fig = plt.scatter(x_50, y_50, c=c, s=100, cmap=plt.cm.get_cmap('RdYlBu_r'))
#cbar = plt.colorbar(fig)
#cbar.set_label('log stellar mass (M$_{\odot}$)', size=20)
#cbar.ax.tick_params(labelsize=20)
plt.vlines(x_50, y_16, y_84, color='gray')
plt.hlines(y_50, x_16, x_84, color='gray')

x = np.arange(0, 40, .1)
plt.plot(x, np.sqrt(1e-3/3*400*x*3), color='k', linestyle='--',
         label='scale height = 400pc, SF duration = 3Gyr')
plt.plot(x, np.sqrt(1e-3/3*180*x*2), color='k', linestyle='-.',
         label='scale height = 180pc, SF duration = 2Gyr')
plt.plot(x, np.sqrt(1e-3/3*80*x*1), color='k', linestyle=':',
         label='scale height = 80pc,   SF duration = 1Gyr')
plt.legend(loc='upper left')


plt.xlim(0, 40)
plt.ylim(0.1, 10)

plt.xlabel('CO velocity dispersion (km/s)', fontsize=20)
plt.ylabel('Correlation length (kpc)', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.savefig(config.savefig_path + 'vel_disp.pdf')
#plt.show()







