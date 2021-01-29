
import config
import constant
from utils import read_CALIFA_catalog
import pandas as pd
import numpy as np


def read_output(name, full_chain):
    array = np.array((full_chain[full_chain[0] == name].values)[0, 1:])
    data_cube = np.reshape(array, [1, config.n_sample, config.n_walker, 4])
    L = np.reshape(data_cube[0, :, :, 2], [1, -1])
    return (np.sqrt(np.percentile(L, 50, axis=1)[0]),
            np.sqrt(np.percentile(L, 16, axis=1)[0]),
            np.sqrt(np.percentile(L, 84, axis=1)[0]))


f = open(config.output_path + '/correlation_length.csv', 'w')
f.write('name,l_50,l_16,l_84\n')

fc = open(config.output_path + '/total_chain_K19N2O2.txt', 'r')

for line in fc.readlines():
    name = line[:line.index(' ')]
    s = line[(1+line.index(' ')):]
    array = np.array([float(n) for n in s.split()])
    data_cube = np.reshape(array, [1, config.n_sample, config.n_walker, 4])
    L = np.reshape(data_cube[0, :, :, 2], [1, -1])
    l_50, l_16, l_84 = (np.sqrt(np.percentile(L, 50, axis=1)[0]),
                        np.sqrt(np.percentile(L, 16, axis=1)[0]),
                        np.sqrt(np.percentile(L, 84, axis=1)[0]))
    f.write('%s,%.3f,%.3f,%.3f\n' %(name, l_50, l_16, l_84))

f.close()
fc.close()


