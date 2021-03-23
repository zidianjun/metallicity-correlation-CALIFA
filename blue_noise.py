
import numpy as np
import matplotlib.pyplot as plt
from cv2 import GaussianBlur
from utils import two_point_correlation
import constant

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
    return (noise.reshape(-1)[mask],
    	    x.reshape(-1)[mask] * kpc_per_pix,
    	    y.reshape(-1)[mask] * kpc_per_pix)

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



