
import constant
import config
import utils

from blue_noise import gen_blue_noise_band
from MCMC import MCMC, KT18_model
import diagnostics

from astropy.io import fits
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import pandas as pd
import time
import uncertainties.unumpy as unp
from uncertainties import ufloat

class EmissionLine(object):
    '''
    A single line with its label, flux, flux uncertainty and rest-frame wavelength.
    '''
    def __init__(self, label, flux, flux_u, unit=1e-16):
        self.label = label
        self.flux = ufloat(flux * unit, flux_u * unit)
        self.rest_wavelength = constant.line_rest_wavelength_dict[label]


class LineTable(object):
    '''
    A table containing all the emission lines.
    '''
    def __init__(self):
        self.label = []
        self.flux = []
        self.rest_wavelength = []

    def read_lines(self, emission_line):
        self.label.append(emission_line.label)
        self.flux.append(emission_line.flux)
        self.rest_wavelength.append(emission_line.rest_wavelength)

    def line_flux(self, label):
        return self.flux[self.label.index(label)]


class Galaxy(object):
    '''
    A class containing galaxy data and properties.
    The data structure is like [[LineTable, ...], [LineTable, ...], ...].
    The eline function serves as a search tool to get the flux of a specified label.
    In morph_par, PA is position angle (degree), b2a is 1 - ellipticity.
    Note that all the data (e.g. eline, x, y0, and met) are 1darray,
        of which the length equals to the number of the pixels of the galaxy.
    They are one-to-one matched and can be easily selected using np.where function (mask).
    Metallicity calculation, two-point correlation function and MCMC model
        are already included in the __init__ function.
    '''
    def __init__(self, gal_name, diag,
                 min_pix=config.min_pix, min_SN=config.min_SN):
        time0 = time.time()
        print "Start analyzing " + gal_name + " with " + diag + " diagnostics."

        self.name = gal_name
        self.min_pix, self.min_SN = min_pix, min_SN
        obj_catalog = utils.read_CALIFA_catalog()
        Pipe3D = fits.open(config.fits_path + gal_name + '.Pipe3D.cube.fits')
        self.eline_data = Pipe3D[3].data
        channel, self.height, self.width = self.eline_data.shape

        FI_data = Pipe3D[1].data[3] # flux intensity
        EBV_data = Pipe3D[1].data[11:13] # dust attenuation
        self.EBV = np.reshape(EBV_data, [2, -1]) / 3.2 # convert from attenuation to E(B-V)
        self.EW = self.eline_data[198].reshape(-1)

        PA, b2a, self.distance, PSF, R_25 = (obj_catalog[obj_catalog['name'] == gal_name].values)[0, 2:7]

        if config.decomposition == 'MA17':
            if gal_name in ['NGC0523', 'UGC04245', 'UGC08107', 'UGC09113', 'NGC7364']:
                pass
            else:
                decom = utils.read_CALIFA_catalog(name='decom.csv')
                ind = (decom.name == gal_name)
                b2a, _, PA = (decom[decom['name'] == gal_name].values)[0, 1:4]

        ic_y, ic_x = int(self.height / 2), int(self.width / 2)
        c_y, c_x = utils.find_max(FI_data, y_range=(ic_y-7, ic_y+7), x_range=(ic_x-8, ic_x+8))

        X, Y = utils.deproject(self.height, self.width, c_y, c_x, PA, b2a, q0=0.)

        self.kpc_per_pix = self.distance * constant.kpc_per_Mpc
        print "1 pixel is %.3fkpc" %(self.kpc_per_pix)
        self.beam = (PSF / np.sqrt(np.log(256)) * self.kpc_per_pix,
                     config.error_PSF / np.sqrt(np.log(256)) * self.kpc_per_pix)

        self.diag = diag
        met, met_u, mask = self.metallicity(self.diag)
        self.mask, self.mask_num = mask, np.sum(mask)
        print "%d pixels left" %(self.mask_num)
        
        if self.mask_num > self.min_pix:

            self.met, self.met_u = met[mask], met_u[mask]
            self.x, self.y = X[mask] * self.kpc_per_pix, Y[mask] * self.kpc_per_pix
            self.rad = np.sqrt(self.x ** 2 + self.y ** 2)

            self.bin_rad, self.bin_met, self.bin_met_u = utils.bootstrap(
                utils.bin_array, (self.met, self.met_u), (self.rad,))
            
            self.fluc, self.fluc_u = utils.step(
                self.rad, self.met, self.met_u,
                self.bin_rad, self.bin_met, self.bin_met_u)

            self.bin_dist, self.bin_ksi, self.bin_ksi_u = utils.bootstrap(
                utils.two_point_correlation, (self.fluc, self.fluc_u), (self.x, self.y))
            
            self.samples, self.par = MCMC(self.bin_dist, self.bin_ksi, self.bin_ksi_u, self.beam)

        else:
            self.met, self.met_u, self.x, self.y = np.nan, np.nan, np.nan, np.nan
            self.rad, self.bin_rad, self.bin_met, self.bin_met_u = np.nan, np.nan, np.nan, np.nan
            self.met_grad, self.met_grad_u = np.nan, np.nan
            self.fluc, self.fluc_u = np.nan, np.nan
            self.bin_dist, self.bin_ksi, self.bin_ksi_u = np.nan, np.nan, np.nan
            self.samples = np.zeros([config.n_sample, config.n_walker, 4])
            self.par = np.ones(4)
        
        time1 = time.time()
        print "Initialization time: %.3fs." %(time1 - time0)
    
    
    def pixel(self, y, x):
        line_table = LineTable()
        for key in config.eline_dict:
            line_table.read_lines(EmissionLine(key,
                                  self.eline_data[config.eline_dict[key]][y][x],
                                  self.eline_data[config.eline_dict[key]+config.diff][y][x]))
        return line_table

    def eline(self, label):
        flux = []
        for y in range(self.height):
            for x in range(self.width):
                flux.append(self.pixel(y, x).line_flux(label))
        return np.array(flux)

    def ratio(self, label_up_list, label_down_list):
        '''
        Dereddened line ratio of two lines using flux and flux uncertainty
        '''
        nan = float('nan')

        label_list = label_up_list + label_down_list
        reddest_wavelength = 0.
        for label in label_list:
            if constant.line_rest_wavelength_dict[label] > reddest_wavelength:
                reddest_wavelength = constant.line_rest_wavelength_dict[label]

        flux_up, flux_down = 0., 0.
        for label in label_up_list:
            uf_array = self.eline(label) 
            uf_array = np.where(uf_array > 0., uf_array, ufloat(nan, nan))
            flux_up += uf_array * diagnostics.relative_dered_factor(self.EBV,
                       constant.line_rest_wavelength_dict[label], reddest_wavelength)
        
        for label in label_down_list:
            uf_array = self.eline(label)
            uf_array = np.where(uf_array > 0., uf_array, ufloat(nan, nan))
            flux_down += uf_array * diagnostics.relative_dered_factor(self.EBV,
                         constant.line_rest_wavelength_dict[label], reddest_wavelength)

        return flux_up / flux_down # uf_array

    def mask_line_flux(self, label_list):
        true_value = True
        for label in label_list:
            uf_array = self.eline(label)
            n, s = diagnostics.split_uf_array(uf_array)
            s = np.where(s > 0., s, -1.)
            true_value = true_value & ((n / s > self.min_SN) & (n > 0.))
        return true_value

    def mask_AGN(self, c=config.AGN_criterion):
        N = unp.log10(self.ratio(['NII6584'], ['Halpha']))
        O = unp.log10(self.ratio(['OIII5007'], ['Hbeta']))
        if c == 'Kewley':
            return ((O < 0.61 / (N - 0.47) + 1.19) & (N < .4))
        elif c == 'Kauffmann':
            return ((O < 0.61 / (N - 0.05) + 1.30))
        else:
            print "AGN criterion must be either 'Kewley' or 'Kauffmann'!"
            return False

    def mask_EW(self, c=config.EW_criterion):
        return np.where(~np.isnan(self.EW), self.EW, 0) < c
    
    def metallicity(self, diag):
        return getattr(diagnostics, diag)(self)
    
    def corr_scale(self, thres=.3):
        if self.bin_ksi is np.nan:
            return np.nan
        for i in range(len(self.bin_ksi)):
            if self.bin_ksi[i+1] <= thres and self.bin_ksi[i] > thres:
                return self.bin_dist[i] + ((self.bin_dist[i+1] - self.bin_dist[i]) /
                       (self.bin_ksi[i+1] - self.bin_ksi[i])) * (thres - self.bin_ksi[i])




class GalaxyFigure(object):
    def __init__(self, galaxy, savefig_path=config.savefig_path):
        self.cmap = plt.cm.get_cmap('RdYlBu_r')
        self.g = galaxy
        '''
        Directly inheriting class Galaxy will double the time
        consumed to fit mcmc.
        '''
        self.savefig_path = savefig_path

    def met_map(self):
        plt.subplots(figsize=(12, 9))
        fig = plt.scatter(self.g.x, self.g.y, c=self.g.met, vmin=8.4, vmax=9.2,
                          s=75, alpha=0.75, edgecolors='face', marker='o', cmap=self.cmap)
        cbar = plt.colorbar(fig)
        cbar.set_label('Metallicity', size=30)
        cbar.ax.tick_params(labelsize=30)
        plt.xlim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.ylim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.xlabel('x (kpc)', fontsize=30)
        plt.ylabel('y (kpc)', fontsize=30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_met_map.pdf')

        plt.subplots(figsize=(12, 9))
        fig = plt.scatter(self.g.x, self.g.y, c=self.g.met_u,
                          marker='o', s=75, alpha=0.75, edgecolors='face', cmap=self.cmap)
        cbar = plt.colorbar(fig)
        cbar.set_label('metallicity uncertainty', size=15)
        cbar.ax.tick_params(labelsize=15)
        plt.xlim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.ylim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.xlabel('x (kpc)', fontsize=15)
        plt.ylabel('y (kpc)', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_met_u_map.pdf')

    def met_fluc(self):
        fig, axes = plt.subplots(2, figsize=(24, 12), sharex=True)
        ax = axes[0]
        ax.errorbar(self.g.rad, self.g.met, yerr=self.g.met_u,
                     marker='.', linestyle='None', color='gray')
        ax.errorbar(self.g.bin_rad, self.g.bin_met, yerr=self.g.bin_met_u,
                     marker='o', linestyle='None', color='k')
        ax.tick_params(axis='both', labelsize=30)
        ax.set_ylabel('Metallicity', fontsize=30)
        ax.set_ylim(8.4, 9.2)

        ax = axes[1]
        ax.errorbar(self.g.rad, self.g.fluc, yerr=self.g.met_u,
                     marker='.', linestyle='None', color='gray')
        ax.axhline(y=0., linestyle='--', color='k')
        ax.tick_params(axis='both', labelsize=30)
        ax.set_ylabel('Metallicity fluctuation', fontsize=30)
        ax.set_xlabel('Radius (kpc)', fontsize=30)
        ax.set_ylim(-.5, .5)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_met_fluc.pdf')

        plt.subplots(figsize=(12, 9))
        max_range = max(abs(min(self.g.fluc)), abs(max(self.g.fluc))) * 1.05
        fig = plt.scatter(self.g.x, self.g.y, c=self.g.fluc, vmin=-.5, vmax=.5,
                          s=75, alpha=0.75, edgecolors='face', marker='o', cmap=self.cmap)
        cbar = plt.colorbar(fig)
        cbar.set_label('Metallicity', size=30)
        cbar.ax.tick_params(labelsize=30)
        plt.xlim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.ylim(-50 * self.g.kpc_per_pix, 50 * self.g.kpc_per_pix)
        plt.xlabel('x (kpc)', fontsize=30)
        plt.ylabel('y (kpc)', fontsize=30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_met_fluc_map.pdf')

    def met_fluc_corr(self):
        plt.subplots(figsize=(8, 4))
        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8)
        
        # For comparison with blue noise.
        if self.g.mask_num > config.min_pix:
            dist, bin_ksi, bin_ksi_u = gen_blue_noise_band(
                self.g.mask, self.g.kpc_per_pix, self.g.beam,
                height=self.g.height, width=self.g.width)
            plt.fill_between(dist, bin_ksi - bin_ksi_u, bin_ksi + bin_ksi_u,
                color='b', alpha=.3, edgecolors='none')
        #print np.sqrt(2 * utils.fit_sigma(dist, bin_ksi)[0][0] / self.g.kpc_per_pix) * 2.354
        
        plt.errorbar(self.g.bin_dist, self.g.bin_ksi, yerr=self.g.bin_ksi_u,
                     linestyle='none', marker='o', color='k', label='galaxy')

        plt.xlim(-.1, 3.6)
        plt.ylim(-.1, 1.1)
        plt.xlabel("$r$ (kpc)", fontsize=15)
        plt.ylabel("$\\xi$($r$)", fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.annotate(self.g.name, xy=(2.7, .8), xytext=(2.7, .8), fontsize=15)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_2p_corr.pdf')

    def mcmc_plot(self):
        plt.subplots(figsize=(8, 4))
        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.8)

        plt.errorbar(self.g.bin_dist, self.g.bin_ksi, yerr=self.g.bin_ksi_u,
                     linestyle='none', marker='o', color='k', label='galaxy')
        
        x = np.arange(.2, 4., .2)
        for i in range(30):
            rs = np.random.randint(config.n_sample)
            rw = np.random.randint(config.n_walker)
            par = self.g.samples[rs, rw]
            plt.plot(np.append([0], x), np.append([1], KT18_model(x, *par[:3]) / par[3]),
                     color='orange', alpha=.3)

        plt.xlim(-.1, 3.6)
        plt.ylim(-.1, 1.1)
        plt.xlabel("$r$ (kpc)", fontsize=15)
        plt.ylabel("$\\xi$($r$)", fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.annotate(self.g.name, xy=(2.7, .8), xytext=(2.7, .8), fontsize=15)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/' + self.g.name + '_' + self.g.diag + '_curves.pdf')

    def example(self):
        plt.figure(figsize=(22, 10))
        plt.subplots_adjust(left=.07, bottom=.1, right=.97, top=.95,
                            wspace=.3, hspace=.0)

        ax = plt.subplot2grid((2, 3), (0, 0))
        fig = ax.scatter(self.g.x, self.g.y, c=self.g.met, vmin=8.4, vmax=9.2,
                         s=10, alpha=0.75, edgecolors='face', marker='o', cmap=self.cmap)
        cbar = plt.colorbar(fig)
        cbar.set_ticks(np.arange(8.4, 9.3, .1))
        cbar.ax.tick_params(labelsize=25)
        ax.set_xticks([])
        ax.set_yticks(np.arange(-10, 15, 5))
        ax.set_ylabel('y (kpc)', fontsize=25)
        ax.tick_params(axis='y', labelsize=25)

        ax = plt.subplot2grid((2, 3), (0, 1), colspan=2)
        ax.errorbar(self.g.rad, self.g.met, yerr=self.g.met_u,
                     marker='.', linestyle='None', color='gray', label='pixels')
        ax.errorbar(self.g.bin_rad, self.g.bin_met, yerr=self.g.bin_met_u,
                     marker='o', linestyle='None', color='k', label='annulus averages')
        ax.legend(loc='lower left')
        ax.set_xticks([])
        ax.tick_params(axis='y', labelsize=25)
        ax.set_ylabel('Metallicity', fontsize=25, labelpad=25)
        ax.set_xlim(-.2, 14)
        ax.set_ylim(8.4, 9.2)

        ax = plt.subplot2grid((2, 3), (1, 0))
        fig = ax.scatter(self.g.x, self.g.y, c=self.g.fluc, vmin=-.5, vmax=.5,
                         s=10, alpha=0.75, edgecolors='face', marker='o', cmap=self.cmap)
        cbar = plt.colorbar(fig)
        cbar.set_ticks(np.arange(-.4, .6, .2))
        cbar.ax.tick_params(labelsize=25)
        ax.set_xticks(np.arange(-10, 15, 5))
        ax.set_yticks(np.arange(-10, 15, 5))
        ax.set_xlabel('x (kpc)', fontsize=25)
        ax.set_ylabel('y (kpc)', fontsize=25)
        ax.tick_params(axis='both', labelsize=25)

        ax = plt.subplot2grid((2, 3), (1, 1), colspan=2)
        ax.errorbar(self.g.rad, self.g.fluc, yerr=self.g.met_u,
                     marker='.', linestyle='None', color='gray')
        ax.axhline(y=0., linestyle='--', color='k')
        ax.tick_params(axis='both', labelsize=25)
        ax.set_xlabel('Radius (kpc)', fontsize=25)
        ax.set_ylabel('Metallicity fluctuation', fontsize=25)
        ax.set_xlim(-.2, 14)
        ax.set_ylim(-.5, .5)

        if self.savefig_path is not None:
            plt.savefig(self.savefig_path + '/fluc_new.pdf')




def analyze(gal_name):

    for diag in config.diag_list:
        galaxy = Galaxy(gal_name, diag)
        samples = galaxy.samples.reshape(-1)
        
        galaxy_figure = GalaxyFigure(galaxy, savefig_path=config.savefigs_path) #
        if samples[0] > 0:
            galaxy_figure.met_map()
            galaxy_figure.met_fluc_corr()
        
        suffix = ''
        if config.AGN_criterion == 'Kauffmann':
            suffix = '_Ka03'
        elif config.adp_bin:
            suffix == '_adp'
        elif config.decomposition == 'MA17':
            suffix = '_MA17'
        
        f = open(config.output_path + '/output/total_chain_' +
                 gal_name + suffix + '.txt', 'a+')
        for i in range(len(samples)):
            if i == 0:
                f.write("%.3f" %(samples[i]))
            else:
                f.write(" %.3f" %(samples[i]))
        f.write("\n")
        f.close()
        

def thres_func(par_tup):
    min_SN, min_pix, name = par_tup
    for i in range(10):
        par = Galaxy(name, 'K19N2O2', min_SN=min_SN, min_pix=min_pix).par
        f = open(config.output_path + name + '_thres.csv', 'a+')
        f.write('%s,%.3f,%.3f,%.3f\n' %(name, min_SN, min_pix, par[2]))
        f.close()

def plot(gal_name):
    galaxy = Galaxy(gal_name, 'K19N2O2')
    galaxy_figure = GalaxyFigure(galaxy, savefig_path=config.savefig_path+'/examples/') #
    galaxy_figure.met_fluc_corr()
    galaxy_figure.mcmc_plot()

def write_corr_scale(gal_name):
    galaxy = Galaxy(gal_name, 'K19N2O2')
    f = open(config.output_path + 'corr_scale.csv', 'a+')
    f.write('%s,%.0f,%.0f\n' %(gal_name, 1e3*galaxy.corr_scale(thres=.5),
                                         1e3*galaxy.corr_scale(thres=.3)))
    f.close()

def write_met_grad(gal_name):
    galaxy = Galaxy(gal_name, 'K19N2O2')
    f = open(config.output_path + 'met_grad.csv', 'a+')
    f.write('%s,%.3f\n' %(gal_name, galaxy.met_grad))
    f.close()


