
from galaxy_metallicity import analyze, thres_func
from utils import read_CALIFA_catalog
import config

import multiprocessing

if __name__ == '__main__':
    
    cs = open(config.output_path + 'corr_scale.csv', 'a+')
    cs.write("ID,name,HW50M,HW30M\n")
    cs.close()
    
    name_list = read_CALIFA_catalog().name
    #name_list = ['NGC0257', 'NGC0776', 'NGC0873', 'NGC1659']  # examples
    cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=cores)
    pool.map(analyze, name_list)
    
    cs = read_CALIFA_catalog(name='corr_scale.csv', path=config.output_path)
    cs.sort_values(by='ID').to_csv(config.output_path + 'corr_scale.csv')


    '''
    tuple_list = []
    for min_SN in [.5, 1., 1.5, 2., 2.5, 3., 3.5]:
        for min_pix in [400, 500, 600, 700]:
            for name in ['NGC0873', 'NGC0237']:
                tuple_list.append((min_SN, min_pix, name))
    pool.map(thres_func, tuple_list)
    '''


