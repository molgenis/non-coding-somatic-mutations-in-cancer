import sys
# import multiprocessing as mp
import pandas as pd
import numpy as np
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from multiprocessing import Pool
import multiprocessing as mp

import ast
import get_data as get_data
import tests as tests


def prep_file(path_file, path_save, type_data, select_chrom):
    df = pd.read_csv(path_file, sep='\t')
    df = df[df['chr'] == select_chrom]

    print(df)
    breast_count_list = list()
    nonbreast_count_list = list()
    for index, row in df.iterrows():
        # print(index)
        if row['cancer_count'] != '-':
            dict_cancer_count = ast.literal_eval(row['cancer_count'])
            total_count = sum(dict_cancer_count.values())
            if 'Breast' in dict_cancer_count:
                breast_count_list.append(dict_cancer_count['Breast'])
                nonbreast_count = total_count - dict_cancer_count['Breast']
                nonbreast_count_list.append(nonbreast_count)
            else:
                breast_count_list.append(0)
                nonbreast_count_list.append(total_count)
        else:
            breast_count_list.append(0)
            nonbreast_count_list.append(0)
            
    df_b_nb = df.iloc[:, 1:5]
    df_b_nb['counts_breast'] = breast_count_list
    df_b_nb['counts_nonbreast'] = nonbreast_count_list
    print(df_b_nb)
    df_b_nb.to_csv(f"{path_save}b_nb_{type_data}_{select_chrom}.tsv", sep='\t', encoding='utf-8', index=False)
    return df_b_nb

def run_all(type_data, path_db, path_save, select_chrom):
    filter_par = False
   
    path_file = f"{path_save}{type_data}_chrALL_num_snps.tsv"

    df_b_nb = prep_file(path_file, path_save, type_data, select_chrom)
    print(df_b_nb)
    print(df_b_nb.head())
    print(set(df_b_nb['counts_breast']))
    print(set(df_b_nb['counts_nonbreast']))    
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par, path_save, path_db)
    # tests_df_all = tests.all_test(df_b_nb, all_num_donor_b, all_num_donor_nb, 'NonCoding', f'{type_data}', path_save, select_chrom)
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, nc_snps_b, nc_snps_nb = get_data.get_noncoding_data(filter_par, path_file, path_db)
    # tests_df_NC = tests.all_test(df_b_nb, nc_snps_b, nc_snps_nb, 'NonCoding_NC', f'{type_data}', path_save, select_chrom)

    parts_df_b_nb = np.array_split(df_b_nb, 20)
    cpus = mp.cpu_count()

    arg_multi_list = []
    for i, df_part in enumerate(parts_df_b_nb):
        arg_multi_list.append((df_part, all_snps_b, all_snps_nb, 'NonCoding', f'{type_data}', path_save, select_chrom, i))
        arg_multi_list.append((df_part, nc_snps_b, nc_snps_nb, 'NonCoding_NC', f'{type_data}', path_save, select_chrom, i))


    pool = Pool(processes=cpus)
    pool.starmap(func=tests.all_test, iterable=arg_multi_list)
    pool.close()
    pool.join()




def main():
    config = get_config()
    path_db = '' #'D:/Hanze_Groningen/STAGE/lastdb/db_laatste_copy.db' #config['database']
    path_save = config['analyse'] #config['analyse'] #'D:/Hanze_Groningen/STAGE/UMAP/'
    type_data = 'UCNE'
    select_chrom = sys.argv[1].replace('chr', '')
    run_all(type_data, path_db, path_save, select_chrom)



if __name__ == '__main__':
    main()

