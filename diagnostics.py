
import constant
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.interpolate import interp1d

# ========== utils for ufloat arrays ========== 

def split_uf_array(uf_array):
    n, s = np.transpose(np.array(map(lambda x: [x.n, x.s], uf_array)))
    return n, s

def merge_uf_array(n_array, s_array):
    assert len(n_array) == len(s_array)
    return np.array(map(lambda n, s: ufloat(n, s), n_array, s_array))

# ========== utils for derendening ==========

def kappa(x, i):
    '''
    From Clayton Cardelli Mathis 1989 law
    '''
    Rv = 3.1  # Clayton Cardelli Mathis 1989
    if i == 1:
        a = 0.574 * x ** 1.61
        b = -0.527 * x ** 1.61
    elif i == 2:
        y = x - 1.82
        a = (1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 +
             0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7)
        b = (1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 -
             0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7)
    elif i == 3:
        a = 1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341)
        b = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263)
    elif i == 4:
        a = (1.752 - 0.316 * x - 0.104 / ((x - 4.67) ** 2 + 0.341) -
             0.04473 * (x - 5.9) ** 2 - 0.009779 * (x - 5.9) ** 3)
        b = (-3.090 + 1.825 * x + 1.206 / ((x - 4.62) ** 2 + 0.263) +
             0.2130 * (x - 5.9) ** 2 - 0.1207 * (x - 5.9) ** 3)
    elif i == 5:
        a = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) ** 2 - 0.070 * (x - 8) ** 3
        b = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) ** 2 + 0.374 * (x - 8) ** 3
    return a * Rv + b

def getfullkappa(wavelength):
    '''
    calculate kappa for full wavelength range
    '''
    x = 1e4 / wavelength # unit of wave should be micron
    k = 0.
    k += kappa(x, 1) * ((x >= 0.3) & (x <= 1.1))
    k += kappa(x, 2) * ((x > 1.1) & (x <= 3.3))
    k += kappa(x, 3) * ((x > 3.3) & (x < 5.9))
    k += kappa(x, 4) * ((x >= 5.9) & (x <= 8.))
    k += kappa(x, 5) * ((x >= 8.) & (x <= 10.))
    return k

def relative_dered_factor(EBV, wavelength, reddest_wavelength):
    EBV = merge_uf_array(EBV[0], EBV[1])
    kappa_a = getfullkappa(wavelength)
    kappa_b = getfullkappa(reddest_wavelength)
    factor = 10 ** (0.4 * (kappa_a - kappa_b) * EBV)
    return factor

# ========== metallicity diagnostics ==========

def PPN2(galaxy):
    N2 = unp.log10(galaxy.ratio(['NII6584'], ['Halpha']))
    met, met_u = split_uf_array(9.37 + 2.03 * N2 + 1.26 * N2 ** 2 + 0.32 * N2 ** 3)
    mask = (galaxy.mask_line_flux(['NII6584', 'Halpha']) &
            galaxy.mask_AGN() & galaxy.mask_EW())
    return (met, np.sqrt(met_u ** 2 + constant.intr_PPN2_error ** 2), mask)

def PPO3N2(galaxy):
    O3 = galaxy.ratio(['OIII5007'], ['Hbeta'])
    N2 = galaxy.ratio(['NII6584'], ['Halpha'])
    O3N2 = unp.log10(O3 / N2)
    met, met_u = split_uf_array(8.73 - 0.32 * O3N2)
    mask = (galaxy.mask_line_flux(['NII6584', 'Halpha', 'OIII5007', 'Hbeta']) &
            galaxy.mask_AGN() & galaxy.mask_EW())
    return (met, np.sqrt(met_u ** 2 + constant.intr_PPO3N2_error ** 2), mask)

def inner(x, y=-3., coeff_logOH=constant.coeff_logOH_N2O2):
    xy = [1, x, y, x*y, x**2, y**2, x*y**2, y*x**2, x**3, y**3]
    m = np.inner(coeff_logOH, xy)
    return ufloat(m.n, m.s)

def K19N2O2(galaxy):
    N2O2_array = unp.log10(galaxy.ratio(['NII6584'], ['OII3727']))

    mask = (galaxy.mask_line_flux(['OII3727', 'NII6584']) &
            galaxy.mask_AGN() & galaxy.mask_EW())
    met = np.array([ufloat(0., 0.)]*len(N2O2_array))
    for index in range(len(N2O2_array)):
        if mask[index]:
            met[index] = inner(N2O2_array[index])
    met, met_u = split_uf_array(met)
    mask = mask & (met > 7.63) & (met < 9.23)
    return (met, np.sqrt(met_u ** 2 + constant.intr_N2O2_error ** 2), mask)




