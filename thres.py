
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


plt.subplots(figsize=(10, 8))
plt.scatter(SN_list+0.03, lcorr_mat[0, 0, :, 0], s=200, marker='^', color='r', label='pixel number = 400')
plt.vlines(SN_list+0.03, lcorr_mat[0, 0, :, 1], lcorr_mat[0, 0, :, 2], color='r', alpha=.5, linewidth=2)
plt.scatter(SN_list+0.01, lcorr_mat[0, 1, :, 0], s=200, marker='s', color='orange', label='pixel number = 500')
plt.vlines(SN_list+0.01, lcorr_mat[0, 1, :, 1], lcorr_mat[0, 1, :, 2], color='orange', alpha=.5, linewidth=2)
plt.scatter(SN_list-0.01, lcorr_mat[0, 2, :, 0], s=200, marker='p', color='cyan', label='pixel number = 600')
plt.vlines(SN_list-0.01, lcorr_mat[0, 2, :, 1], lcorr_mat[0, 2, :, 2], color='cyan', alpha=.5, linewidth=2)
plt.scatter(SN_list-0.03, lcorr_mat[0, 3, :, 0], s=200, marker='h', color='b', label='pixel number = 700')
plt.vlines(SN_list-0.03, lcorr_mat[0, 3, :, 1], lcorr_mat[0, 3, :, 2], color='b', alpha=.5, linewidth=2)
plt.xlim(0, max(SN_list) + .5)
plt.ylim(-.5, 6.5)
plt.xlabel('S/N threshold', fontsize=20)
plt.ylabel('Correlation length', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.annotate('NGC0873', xy=(2.5, 3.5), xytext=(2.5, 3.5), fontsize=20)
plt.legend(loc='upper right', fontsize=20)
choice = lcorr_mat[0, 1, 2, :]
plt.scatter(2+0.01, choice[0], s=200, marker='s', color='orange', edgecolor='k')
plt.vlines(2+0.01,choice[1], choice[2], color='k')
plt.savefig(config.savefig_path + 'thres_NGC0873.pdf')

plt.subplots(figsize=(10, 8))
plt.scatter(SN_list+0.03, lcorr_mat[1, 0, :, 0], s=200, marker='^', color='r', label='pixel number = 400')
plt.vlines(SN_list+0.03, lcorr_mat[1, 0, :, 1], lcorr_mat[1, 0, :, 2], color='r', alpha=.5, linewidth=2)
plt.scatter(SN_list+0.01, lcorr_mat[1, 1, :, 0], s=200, marker='s', color='orange', label='pixel number = 500')
plt.vlines(SN_list+0.01, lcorr_mat[1, 1, :, 1], lcorr_mat[1, 1, :, 2], color='orange', alpha=.5, linewidth=2)
plt.scatter(SN_list-0.01, lcorr_mat[1, 2, :, 0], s=200, marker='p', color='cyan', label='pixel number = 600')
plt.vlines(SN_list-0.01, lcorr_mat[1, 2, :, 1], lcorr_mat[1, 2, :, 2], color='cyan', alpha=.5, linewidth=2)
plt.scatter(SN_list-0.03, lcorr_mat[1, 3, :, 0], s=200, marker='h', color='b', label='pixel number = 700')
plt.vlines(SN_list-0.03, lcorr_mat[1, 3, :, 1], lcorr_mat[1, 3, :, 2], color='b', alpha=.5, linewidth=2)
plt.xlim(0, max(SN_list) + .5)
plt.ylim(-.5, 6.5)
plt.xlabel('S/N threshold', fontsize=20)
plt.ylabel('Correlation length', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.annotate('NGC0237', xy=(2.5, 3.5), xytext=(2.5, 3.5), fontsize=20)
plt.legend(loc='upper right', fontsize=20)
choice = lcorr_mat[1, 1, 2, :]
plt.scatter(2+0.01, choice[0], s=200, marker='s', color='orange', edgecolor='k')
plt.vlines(2+0.01, choice[1], choice[2], color='k')
plt.savefig(config.savefig_path + 'thres_NGC0237.pdf')

#plt.show()




