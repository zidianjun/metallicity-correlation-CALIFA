
# data path
proj_path = './'
#proj_path = '/home/liz/project/metallicity/'
obj_path = proj_path + '/data/'
fits_path = proj_path + '/data/fits/'
savefig_path = proj_path + '/figure/'
savefigs_path = proj_path + '/figures/'
output_path = proj_path + '/output/'

# CALIFA columns
eline_dict = {'OII3727': 0, 'OIII5007': 26, 'OIII4959': 27, 'Hbeta': 28,
              'Halpha': 45, 'NII6584': 46, 'NII6549': 47, 'SII6717': 49, 'SII6731': 50}

diff = 204

# metallicity diagnostics
diag_list = ['PPN2', 'PPO3N2', 'K19N2O2']

# adaptive bin width
adp_bin = False

# galaxy properties
min_SN = 2.
min_pix = 500
q0 = .13  # .13 for CALIFA DR3 thick-disk de-projection
EW_criterion = -6  # Sanchez
AGN_criterion = 'Kewley'  # or Kauffmann


# PSF is around 2.5".
error_PSF = .15

# MCMC
n_walker = 100
n_step = 500
n_sample = 150



