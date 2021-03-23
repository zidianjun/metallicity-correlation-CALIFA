
import config
from utils import read_CALIFA_catalog
import numpy as np
import os

obj_catalog = read_CALIFA_catalog()

suffix = ''
if config.AGN_criterion == 'Kauffmann':
    suffix = '_Ka03'
elif config.adp_bin:
    suffix == '_adp'
elif config.decomposition == 'MA17':
    suffix = '_MA17'

PPN2 = open(config.output_path + '/total_chain_PPN2' + suffix + '.txt', 'w')
PPO3N2 = open(config.output_path + '/total_chain_PPO3N2' + suffix + '.txt', 'w')
K19N2O2 = open(config.output_path + '/total_chain_K19N2O2' + suffix + '.txt', 'w')
for name in obj_catalog['name']:
	path = config.output_path + '/output/total_chain_' + name + suffix + '.txt'
	f = open(path, 'r')
	PPN2.write("%s " %(name))
	PPO3N2.write("%s " %(name))
	K19N2O2.write("%s " %(name))
	for i, line in enumerate(f.readlines()):
		if i == 0:
			PPN2.write(line)
		if i == 1:
			PPO3N2.write(line)
		if i == 2:
			K19N2O2.write(line)
	f.close()


PPN2.close()
PPO3N2.close()
K19N2O2.close()
