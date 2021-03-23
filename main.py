
from galaxy_metallicity import analyze, thres_func, plot, write_corr_scale, write_met_grad
from utils import read_CALIFA_catalog

import multiprocessing

if __name__ == '__main__':

    name_list = read_CALIFA_catalog(name='name_list.csv').name
    #name_list = ['NGC0257', 'NGC0776', 'NGC0873', 'NGC1659']  # examples
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    pool.map(analyze, name_list)
    '''
    tuple_list = []
    for min_SN in [.5, 1., 1.5, 2., 2.5, 3., 3.5]:
        for min_pix in [400, 500, 600, 700]:
            for name in ['NGC0873', 'NGC0237']:
                tuple_list.append((min_SN, min_pix, name))
    pool.map(thres_func, tuple_list)
    '''
    #for name in name_list:
    #    write_met_grad(name)


