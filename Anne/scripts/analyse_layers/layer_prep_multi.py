#!/usr/bin/env python3

# Imports
import sys
import numpy as np

sys.path.append(
    '/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from multiprocessing import Pool
import multiprocessing as mp
from prep_file_for_analyse import prep_file_multi
import get_data as get_data
import tests as tests


def layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, nc_c):
    """
    Ensures that the various tests are run multiprocess
    :param type_data: Type of data/name of the layer.
    :param path_db:  Path to the database
    :param path_save:     Path where the files should be saved
    :param select_chrom:  The selected chromosome
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param nc_c:  Is the data type noncoding, coding or a mix
    :return:    
    """
    if type_data in ['after', 'before']:
        path_file = f"{path_save}ALL_gene_{type_data}_2000_250.tsv"
    else:
        path_file = f"{path_save}{type_data}_chrALL_num_snps.tsv"
    # Call prep_file_multi
    df_b_nb = prep_file_multi(path_file, path_save, type_data, select_chrom)
    # Get all data and noncoding data   
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(
        filter_par, path_save, path_db)
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, nc_snps_b, nc_snps_nb = get_data.get_noncoding_data(
        filter_par, path_file, path_db)
    # Multiprocess
    parts_df_b_nb = np.array_split(df_b_nb, 20)
    cpus = mp.cpu_count()
    arg_multi_list = []
    for i, df_part in enumerate(parts_df_b_nb):
        arg_multi_list.append((df_part, all_snps_b, all_snps_nb, nc_c, f'{type_data}', path_save, select_chrom, i))
        arg_multi_list.append(
            (df_part, nc_snps_b, nc_snps_nb, f'{nc_c}_NC', f'{type_data}', path_save, select_chrom, i))
    pool = Pool(processes=cpus)
    pool.starmap(func=tests.all_test, iterable=arg_multi_list)
    pool.close()
    pool.join()


def main():
    # Call get_config
    config = get_config('calculon')
    # When a database path is listed here, files still need to be created.
    # When this is read, it is assumed that no files need to be created.
    path_db = ''
    path_save = config['analyse']
    select_chrom = sys.argv[1].replace('chr', '')
    filter_par = False
    # Run for different layers
    type_data = 'DNase'
    layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, 'NonCoding')
    type_data = 'TFBS'
    layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, 'NonCoding')
    type_data = 'UCNE'
    layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, 'NonCoding')
    type_data = 'before'
    layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, 'NonCoding_Coding')
    type_data = 'after'
    layers_run_multi(type_data, path_db, path_save, select_chrom, filter_par, 'NonCoding_Coding')


if __name__ == '__main__':
    main()
