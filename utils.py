
import constant
import config

import numpy as np
from pandas import read_csv
from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.special import jv

# ========== read table ========== 

def read_CALIFA_catalog(name='CALIFA_galaxy.csv', path=config.obj_path):
    return read_csv(path + name)

def open_full_chain(diag='K19N2O2', path=config.output_path, suffix=''):
    return open(path + '/total_chain_' + diag + suffix + '.txt', 'r')

# ========== multi-parameter fitting ========== 

def grad(x, y, e):

    def linear(x, a, b):
        return a * x + b

    p0 = [-0.05, 9]

    popt, pcov = curve_fit(linear, x, y, sigma=e)
    s_sq = ((linear(x, *popt) - y) ** 2).sum() / (len(x) - len(p0))
    return popt[0], np.sqrt(np.diag(pcov * s_sq))[0]

# ========== deprojection ========== 

def find_max(matrix, y_range=None, x_range=None):
    height, width = matrix.shape
    x_min, x_max = x_range if x_range is not None else (0, width)
    y_min, y_max = y_range if y_range is not None else (0, height)
    mat = matrix[y_min:y_max, x_min:x_max]
    index = np.where(mat == np.max(mat))
    return index[0][0] + y_min, index[1][0] + x_min

def deproject(height, width, y_c, x_c, PA, b2a, q0=config.q0):
    '''
    Deproject the galaxy coordinates using rotation matrix.
    '''
    cosi = np.sqrt((b2a ** 2 - q0 ** 2) / (1 - q0 ** 2)) if b2a > q0 else 0.
    theta = PA * np.pi / 180
    dep_mat = np.array([[np.cos(theta), np.sin(theta)],
                        [-np.sin(theta) / cosi, np.cos(theta) / cosi]])
    
    x0, y0 = np.meshgrid(range(width), range(height))
    x0, y0 = np.squeeze(np.reshape(x0, [1, -1])), np.squeeze(np.reshape(y0, [1, -1]))
    xy_mat = np.stack([x0 - x_c, y0 - y_c], axis=0)
    X, Y = np.dot(dep_mat, xy_mat)
    return (X, Y)

# ========== utils for 1d np arrays ==========


def bin_array(f, r, bin_size=.2, adp=config.adp_bin): # phase=0
    if adp:
        hist, bin_edge = np.histogram(r, bins='auto')
    else:
        bin_edge = np.arange(0, max(r) + bin_size, bin_size)
    stat = binned_statistic(r, f, bins=bin_edge)
    mask = ~np.isnan(stat.statistic)
    bin_r = stat.bin_edges[:-1][mask]
    bin_f = stat.statistic[mask]
    return bin_r, bin_f

def step(rad, met, met_u, bin_rad, bin_met, bin_met_u):
    rad_matrix = abs(np.subtract.outer(rad, bin_rad))
    min_value = np.expand_dims(np.min(rad_matrix, axis=1), axis=1)
    x, y = np.where(rad_matrix == min_value) # np.where will return ALL the min values!
    y = y[np.insert(np.diff(x) != 0, 0, True)] # If having recurring items, this removes them.
    step_func = bin_met[y]
    fluc = met - step_func
    fluc_u = np.sqrt(met_u ** 2 + bin_met_u[y] ** 2)
    return fluc, fluc_u

def two_point_correlation(f, x, y, bin_size=.2, short=20):
    SCORR = np.outer(f, f)
    DX, DY = np.subtract.outer(x, x), np.subtract.outer(y, y)
    DIST = np.sqrt(DX ** 2 + DY ** 2)

    scorr = np.squeeze(np.reshape(SCORR, [1, -1]))
    dist = np.squeeze(np.reshape(DIST, [1, -1]))
    mean2, sigma2 = np.mean(f) ** 2, np.std(f) ** 2 # mean is 0
    
    stat = binned_statistic(dist, scorr, bins=np.arange(0, max(dist) + bin_size, bin_size))
    mask = ~np.isnan(stat.statistic)
    bin_d = stat.bin_edges[:-1][mask][:short]
    bin_s = (stat.statistic[mask][:short] - mean2) / sigma2

    return bin_d, bin_s

def bootstrap(func, bootstrap_args, fixed_args, times=20):
    value_b = None
    for time in range(times):
        bootstrap_value = np.random.normal(*bootstrap_args)
        x, value = func(bootstrap_value, *fixed_args)
        value = np.expand_dims(value, axis=0)
        value_b = np.concatenate((value_b, value), axis=0) if value_b is not None else value
        #print "Bootstrap for %s #%d finished." %(func.__name__, time)
    value_b_mean = np.mean(value_b, axis=0)
    value_b_std = np.std(value_b, axis=0)
    return x, value_b_mean, value_b_std


