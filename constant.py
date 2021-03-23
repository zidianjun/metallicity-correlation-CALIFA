
import numpy as np

# In CALIFA survey, 1 pixel equals to 1 arcsec, thus, 4.848e-6, 4.848e-3 kpc/Mpc
kpc_per_Mpc = 4.848e-3

line_rest_wavelength_dict = {
	'OII3727': 3727.092,
	'Hbeta': 4862.683,
	'OIII4959': 4960.295,
	'OIII5007': 5008.240,
	'NII6549': 6549.840,
	'Halpha': 6564.610,
	'NII6584': 6585.230,
	'SII6717': 6718.294,
	'SII6731': 6732.674
}


coeff_logOH_N2O2 = [9.4774, 1.1789, 0.5092, 0.6867, 0.2814, 0.1617, 0.1184, 0.1202, 0.2292, 0.0165]

intr_PPN2_error = 0.021
intr_PPO3N2_error = 0.012
intr_N2O2_error = np.log10(1 + 0.0265)