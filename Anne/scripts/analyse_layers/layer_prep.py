#!/usr/bin/env python3

# Imports
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from prep_file_for_analyse import prep_file
import get_data as get_data
import tests as tests


def layers_run(type_data, path_db, filter_par, path_save):
    """

    :param : 
    :param :  
    :param :        
    :return:    
    """
    if type_data in ['after', 'before']:
        path_file = f"{path_save}ALL_gene_{type_data}_2000_250.tsv"
        df_b_nb = prep_file(path_file, path_save, type_data)    
        # Call run_all for before and after
        all_nc = 'NonCoding_Coding'
        all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par, path_save, path_db)
        tests_df_all = tests.all_test(df_b_nb, all_num_donor_b, all_num_donor_nb, all_nc, f'{type_data}', path_save, 'chr0', '0')
    else:
        path_file = f"{path_save}{type_data}_chrALL_num_snps.tsv"
        df_b_nb = prep_file(path_file, path_save, type_data) 
        all_nc = 'NonCoding2'
        all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par, path_save, path_db)
        tests_df_all = tests.all_test(df_b_nb, all_num_donor_b, all_num_donor_nb, all_nc, f'{type_data}', path_save, 'chr0', '0')
        all_nc = 'NonCoding_NC2'
        noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(filter_par, path_save, path_db)
        tests_df_NC = tests.all_test(df_b_nb, all_snps_b, all_snps_nb, all_nc, f'{type_data}', path_save, 'chr0', '0')    
    


def main():
    # Call get_config
    config = get_config('gearshift')
    # When a database path is listed here, files still need to be created.
    # When this is read, it is assumed that no files need to be created.
    path_db = '' 
    filter_par = False
    path_save = config['analyse']
    type_data = 'DNase'
    layers_run(type_data, path_db, filter_par, path_save)
    type_data = 'TFBS'
    layers_run(type_data, path_db, filter_par, path_save)
    type_data = 'UCNE'
    layers_run(type_data, path_db, filter_par, path_save)
    type_data = 'before'
    layers_run(type_data, path_db, filter_par, path_save)
    type_data = 'after'
    layers_run(type_data, path_db, filter_par, path_save)




if __name__ == '__main__':
    main()

