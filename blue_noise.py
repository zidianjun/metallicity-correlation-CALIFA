
import numpy as np
import matplotlib.pyplot as plt
from cv2 import GaussianBlur
from utils import two_point_correlation
import constant
import matplotlib.pyplot as plt

def gen_blue_noise(mask, kpc_per_pix, beam, height=73, width=78, fsc=5):
    '''
    Height and width are in the unit of pixel.
    1 pixel in CALIFA is 1 arcsec.
    The fine-structure constant(fsc) is serving to produce a higher-resolution sky map.
    The kernel size must be odd and here the default value is 15.
    '''
    noise = GaussianBlur(np.random.normal(0, 1, (width * fsc, height * fsc)),
                         (15, 15), beam[0] / kpc_per_pix * fsc)[::fsc, ::fsc]
    y, x = np.meshgrid(np.arange(height), np.arange(width))
    return (np.squeeze(np.reshape(noise, [1, -1]))[mask],
    	    np.squeeze(np.reshape(x, [1, -1]))[mask] * kpc_per_pix,
    	    np.squeeze(np.reshape(y, [1, -1]))[mask] * kpc_per_pix)

def gen_blue_noise_band(mask, kpc_per_pix, beam, height=73, width=78, times=10):
    ksi_b = None
    for time in range(times):
        x, ksi = two_point_correlation(*gen_blue_noise(mask, kpc_per_pix, beam,
                                                       height=height, width=width))
        ksi = np.expand_dims(ksi, axis=0)
        ksi_b = np.concatenate((ksi_b, ksi), axis=0) if ksi_b is not None else ksi
        print "Bootstrap for blue noise #%d finished." %(time)
    ksi_b_mean = np.mean(ksi_b, axis=0)
    ksi_b_std = np.std(ksi_b, axis=0)
    return x, ksi_b_mean, ksi_b_std



'''

def gen_noise(scale=1, height=50, width=50, beam=3, fsc=5):
    noise = GaussianBlur(np.random.normal(0, 1, (width * fsc, height * fsc)),
                         (15, 15), beam * fsc)[::fsc, ::fsc]
    y, x = np.meshgrid(np.arange(height), np.arange(width))
    return (np.squeeze(np.reshape(noise, [1, -1])),
            np.squeeze(np.reshape(x, [1, -1])) * scale,
            np.squeeze(np.reshape(y, [1, -1])) * scale)

def gen_noise_band(beam, scale=1, short=30, times=10):
    ksi_b = None
    for time in range(times):
        x, ksi = two_point_correlation(*gen_noise(beam=beam, scale=scale), short=short)
        ksi = np.expand_dims(ksi, axis=0)
        ksi_b = np.concatenate((ksi_b, ksi), axis=0) if ksi_b is not None else ksi
        print "Bootstrap for blue noise #%d finished." %(time)
    ksi_b_mean = np.mean(ksi_b, axis=0)
    ksi_b_std = np.std(ksi_b, axis=0)
    return x, ksi_b_mean, ksi_b_std

if __name__ == '__main__':
    temp = np.random.normal(0, 1, (50, 50))
    plt.figure()
    plt.imshow(GaussianBlur(temp, (15, 15), 1e-6))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.figure()
    plt.imshow(GaussianBlur(temp, (15, 15), 1.))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.figure()
    dist, bin_ksi, bin_ksi_u = gen_noise_band(1e-6, scale=.1, short=-1)
    plt.fill_between(dist, bin_ksi - bin_ksi_u, bin_ksi + bin_ksi_u,
                         color='b', alpha=.3, edgecolors='none')
    plt.xlabel("$x$", fontsize=20)
    plt.ylabel("$\\xi$($x$)", fontsize=20)
    plt.suptitle('two-point correlation', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(-.1, 5.1)
    plt.ylim(-.1, 1.1)
    plt.figure()
    dist, bin_ksi, bin_ksi_u = gen_noise_band(1.)
    plt.fill_between(dist, bin_ksi - bin_ksi_u, bin_ksi + bin_ksi_u,
                         color='b', alpha=.3, edgecolors='none')
    plt.xlabel("$x$", fontsize=20)
    plt.ylabel("$\\xi$($x$)", fontsize=20)
    plt.suptitle('two-point correlation', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(-.1, 5.1)
    plt.ylim(-.1, 1.1)
    plt.show()
'''



