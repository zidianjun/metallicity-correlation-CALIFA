
from galaxy_metallicity import analyze, thres_func, plot, write_corr_scale
from utils import read_CALIFA_catalog

import multiprocessing

if __name__ == '__main__':

    name_list = read_CALIFA_catalog().name
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    pool.map(analyze, name_list)


