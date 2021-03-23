
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from utils import read_CALIFA_catalog
import config


thres0237 = read_CALIFA_catalog(name='NGC0237_thres.csv',
            path=config.output_path+'/thres/')
thres0873 = read_CALIFA_catalog(name='NGC0873_thres.csv',
            path=config.output_path+'/thres/')

SN_list = np.array([3., 2.5, 2., 1.5, 1., .5])
pix_list = [400, 500, 600, 700]

lcorr_mat = np.zeros((2, len(pix_list), len(SN_list), 3))

for i, pix in enumerate(pix_list):
    for j, SN in enumerate(SN_list):
        temp = thres0873[(thres0873.pix == pix) & (thres0873.SN == SN)].lcorr
        x, dx = np.median(temp), np.std(temp) / np.sqrt(len(temp))
        median, ste = np.sqrt(x), 1./2/np.sqrt(x) * dx
        lcorr_mat[0, i, j, 0] = median
        lcorr_mat[0, i, j, 1] = max(0, median - ste)
        lcorr_mat[0, i, j, 2] = median + ste
        temp = thres0237[(thres0237.pix == pix) & (thres0237.SN == SN)].lcorr
        x, dx = np.median(temp), np.std(temp) / np.sqrt(len(temp))
        median, ste = np.sqrt(x), 1./2/np.sqrt(x) * dx
        lcorr_mat[1, i, j, 0] = median
        lcorr_mat[1, i, j, 1] = max(0, median - ste)
        lcorr_mat[1, i, j, 2] = median + ste

plt.figure(figsize=(8, 12))
plt.subplots_adjust(left=.10, bottom=.10, right=.95, top=.95, hspace=.0)

ax = plt.subplot(211)
ax.scatter(SN_list+0.03, lcorr_mat[0, 0, :, 0], s=200, marker='^', color='r', label='pixel number = 400')
ax.vlines(SN_list+0.03, lcorr_mat[0, 0, :, 1], lcorr_mat[0, 0, :, 2], color='r', alpha=.5, linewidth=2)
ax.scatter(SN_list+0.01, lcorr_mat[0, 1, :, 0], s=200, marker='s', color='orange', label='pixel number = 500')
ax.vlines(SN_list+0.01, lcorr_mat[0, 1, :, 1], lcorr_mat[0, 1, :, 2], color='orange', alpha=.5, linewidth=2)
ax.scatter(SN_list-0.01, lcorr_mat[0, 2, :, 0], s=200, marker='p', color='cyan', label='pixel number = 600')
ax.vlines(SN_list-0.01, lcorr_mat[0, 2, :, 1], lcorr_mat[0, 2, :, 2], color='cyan', alpha=.5, linewidth=2)
ax.scatter(SN_list-0.03, lcorr_mat[0, 3, :, 0], s=200, marker='h', color='b', label='pixel number = 700')
ax.vlines(SN_list-0.03, lcorr_mat[0, 3, :, 1], lcorr_mat[0, 3, :, 2], color='b', alpha=.5, linewidth=2)
ax.set_xlim(0, max(SN_list) + .5)
ax.set_ylim(-.5, 6.5)
ax.set_xlabel('S/N threshold', fontsize=20)
ax.set_ylabel('Correlation length', fontsize=20)
ax.tick_params(axis='both', labelsize=20, labelbottom=False)
ax.annotate('NGC0873', xy=(2.5, 3.5), xytext=(2.5, 3.5), fontsize=18)
ax.legend(loc='upper right', fontsize=20, prop={'size': 18})
choice = lcorr_mat[0, 1, 2, :]
ax.scatter(2+0.01, choice[0], s=200, marker='s', color='orange', edgecolor='k')
ax.vlines(2+0.01,choice[1], choice[2], color='k')

ax = plt.subplot(212)
ax.scatter(SN_list+0.03, lcorr_mat[1, 0, :, 0], s=200, marker='^', color='r', label='pixel number = 400')
ax.vlines(SN_list+0.03, lcorr_mat[1, 0, :, 1], lcorr_mat[1, 0, :, 2], color='r', alpha=.5, linewidth=2)
ax.scatter(SN_list+0.01, lcorr_mat[1, 1, :, 0], s=200, marker='s', color='orange', label='pixel number = 500')
ax.vlines(SN_list+0.01, lcorr_mat[1, 1, :, 1], lcorr_mat[1, 1, :, 2], color='orange', alpha=.5, linewidth=2)
ax.scatter(SN_list-0.01, lcorr_mat[1, 2, :, 0], s=200, marker='p', color='cyan', label='pixel number = 600')
ax.vlines(SN_list-0.01, lcorr_mat[1, 2, :, 1], lcorr_mat[1, 2, :, 2], color='cyan', alpha=.5, linewidth=2)
ax.scatter(SN_list-0.03, lcorr_mat[1, 3, :, 0], s=200, marker='h', color='b', label='pixel number = 700')
ax.vlines(SN_list-0.03, lcorr_mat[1, 3, :, 1], lcorr_mat[1, 3, :, 2], color='b', alpha=.5, linewidth=2)
ax.set_xlim(0, max(SN_list) + .5)
ax.set_ylim(-.5, 6.5)
ax.set_xlabel('S/N threshold', fontsize=20)
ax.set_ylabel('Correlation length', fontsize=20)
ax.tick_params(axis='both', labelsize=20)
ax.annotate('NGC0237', xy=(2.5, 3.5), xytext=(2.5, 3.5), fontsize=18)
ax.legend(loc='upper right', fontsize=20, prop={'size': 18})
choice = lcorr_mat[1, 1, 2, :]
ax.scatter(2+0.01, choice[0], s=200, marker='s', color='orange', edgecolor='k')
ax.vlines(2+0.01, choice[1], choice[2], color='k')


plt.savefig(config.savefig_path + 'thres.pdf')
#plt.show()




