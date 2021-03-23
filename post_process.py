
import config
import constant
from utils import read_CALIFA_catalog, open_full_chain
import pandas as pd
import numpy as np

class DataHub(object):
    def __init__(self, diag='K19N2O2',
                 basic=read_CALIFA_catalog(),
                 name_list=read_CALIFA_catalog(name='name_list.csv'),
                 PSF=read_CALIFA_catalog(name='PSF.csv'),
                 cs=read_CALIFA_catalog(name='corr_scale.csv',
                                        path=config.output_path)):
        self.diag = diag
        self.basic = basic
        self.mask = self.basic.name.isin(list(name_list.name))
        self.PSF = PSF
        self.cs = cs

    def col(self, col_name):
        if col_name == 'FWHM':
            return np.array(self.PSF.FWHM[self.mask])
        else:
            return np.array(self.basic[col_name][self.mask])

    def corr_len(self):
        suffix = '' if self.diag == 'K19N2O2' else '_' + self.diag
        lcorr = read_CALIFA_catalog(name='correlation_length' +
                                    suffix +'.csv', path=config.output_path)
        return (np.array(lcorr['l_50'][self.mask]),
                np.array(lcorr['l_16'][self.mask]),
                np.array(lcorr['l_84'][self.mask]))

    def corr_scale(self):
        return (np.array(self.cs['HW50M'][self.mask]),
                np.array(self.cs['HW30M'][self.mask]))

    def rand_corr_len(self, diag='K19N2O2', suffix=''):
        array = []
        f = open_full_chain(diag=diag, suffix=suffix)

        for line in f.readlines():
            name = line[:line.index(' ')]
            if name in list(self.basic.name[self.mask]):
                s = line[(1+line.index(' ')):]
                chain = np.array([float(n) for n in s.split()])
                data_cube = np.reshape(chain, [1, config.n_sample, config.n_walker, 4])
                
                L = 0
                while L == 0:
                    rs = np.random.randint(config.n_sample)
                    rw = np.random.randint(config.n_walker)
                    L = np.mean(data_cube[0, rs, rw, 2])
                
                array.append(np.sqrt(L))

        return np.array(array)

        




