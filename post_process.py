
import config
import constant
from utils import read_CALIFA_catalog, open_full_chain
import pandas as pd
import numpy as np

class DataHub(object):
    def __init__(self, diag='K19N2O2',
                 basic=read_CALIFA_catalog(),
                 cs=read_CALIFA_catalog(name='corr_scale.csv',
                                        path=config.output_path)):
        self.diag = diag
        self.basic = basic
        self.cs = cs

    def col(self, col_name):
        return np.array(self.basic[col_name])

    def corr_len(self):
        suffix = '' if self.diag == 'K19N2O2' else '_' + self.diag
        lcorr = read_CALIFA_catalog(name='correlation_length' +
                                    suffix +'.csv', path=config.output_path)
        return (np.array(lcorr['l_50']),
                np.array(lcorr['l_16']),
                np.array(lcorr['l_84']))

    def corr_scale(self):
        return (np.array(self.cs['HW50M']),
                np.array(self.cs['HW30M']))

    def rand_corr_len(self, diag='K19N2O2', suffix=''):
        array = []
        f = open_full_chain(diag=diag, suffix=suffix)

        for line in f.readlines():
            name = line[:line.index(' ')]
            
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

        




