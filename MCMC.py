
import numpy as np
from scipy.integrate import quad
from scipy.special import jv
from scipy.interpolate import interp2d
import config
import emcee
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
#import corner


def KT18(alpha, beta):
    return (2. / np.log(1 + beta / alpha) * quad(lambda x:
            np.exp(-alpha * x**2) * (1 - np.exp(-beta * x**2)) * jv(0, x)/x, 0, np.inf)[0])

def KT18_model(x_array, sigma, x0, KappaTstar):
    res = []
    for x in x_array:
        alpha = (sigma ** 2 / 2 + x0 ** 2) / x ** 2
        beta = 2 * KappaTstar / x ** 2
        if -12 < np.log(alpha) < 6 and -12 < np.log(beta) < 14:
            res.append(KT18(alpha, beta))
        else:
            res.append(np.inf)
    return np.array(res)

def log_likelihood(theta, x, y, yerr):
    sigma, x0, KappaTstar, f = theta
    model = KT18_model(x, sigma, x0, KappaTstar)
    log_prob = -.5 * np.sum((y - model / f) ** 2 / yerr ** 2 + np.log(yerr ** 2))
    return log_prob

def log_prior(theta, beam):
    sigma, x0, KappaTstar, f = theta
    sigma_0, sigma_u = beam
    if x0 > 0 and KappaTstar > 0 and f > 1 and sigma > 0:
        return -.5 * ((sigma - sigma_0) ** 2 / sigma_u ** 2 + np.log(sigma_u ** 2 * f ** 2))
    return -np.inf

def log_probability(theta, x, y, yerr, beam):
    p = log_prior(theta, beam)
    l = log_likelihood(theta, x, y, yerr)
    if not np.isfinite(p) or np.isnan(p) or not np.isfinite(l) or np.isnan(l):
        return -np.inf
    return p + l

def MCMC(x, y, yerr, beam):
    pos = np.array([beam[0], .01, 1., 1.3]) + 1e-4 * np.random.randn(config.n_walker, 4)
    nwalkers, ndim = pos.shape
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x[1:], y[1:], yerr[1:], beam))
    # The first term of y should always be zero and will never be affected by the factor of f.
    sampler.run_mcmc(pos, config.n_step, progress=True)
    samples = sampler.get_chain()
    flat_samples = np.reshape(samples[-config.n_sample:, :, :], [-1, 4])
    labels = ["$\sigma_{\mathrm{beam}}$", "$\sigma_{\mathrm{inj}}$", "$\kappa t_*$", "$f$"]

    '''
    fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    plt.savefig(config.savefig_path + '/chains.pdf')
    
    corner.corner(flat_samples, labels=labels)
    plt.savefig(config.savefig_path + '/corner.pdf')
    '''

    print(np.percentile(flat_samples, 50, axis=0))
    print(np.percentile(flat_samples, 16, axis=0))
    print(np.percentile(flat_samples, 84, axis=0))
    
    return samples[-config.n_sample:, :, :], np.percentile(flat_samples, 50, axis=0)
    

