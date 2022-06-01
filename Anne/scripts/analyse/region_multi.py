import sys
# import multiprocessing as mp
import pandas as pd
import numpy as np
# sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
sys.path.append('/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')

from config import get_config
from multiprocessing import Pool
import multiprocessing as mp

import get_data as get_data
import tests as tests

def make_merge_df(path_breast, path_nonbreast, select_chrom):
    colnames=['index', 'counts_breast', 'chr', 'start_region', 'stop_region']
    breast = pd.read_csv(path_breast, sep='\t', header=None, names=colnames)
    breast = breast[breast['chr'] == select_chrom]
    breast.sort_values(['chr', 'start_region'], inplace=True)
    
    colnames=['index', 'counts_nonbreast', 'chr', 'start_region', 'stop_region']
    nonbreast = pd.read_csv(path_nonbreast, sep='\t', header=None, names=colnames)
    nonbreast = nonbreast[nonbreast['chr'] == select_chrom]
    nonbreast.sort_values(['chr', 'start_region'], inplace=True)
    
    merged_df = breast.merge(nonbreast, on=['chr', 'start_region', 'stop_region'], how='outer')
    merged_df.drop(columns=['index_x', 'index_y'], inplace=True)
    return merged_df


def make_df_2000(merge_df):
    df = pd.DataFrame(columns=['counts_breast', 'chr', 'start_region', 'stop_region', 'counts_nonbreast', 'snps_b_double', 'snps_nb_double'])

    # chrom = ['chr4']
    for i in list(set(merge_df['chr'])):
        select_df = merge_df[merge_df['chr'] == i]
        select_df['snps_b_double'] = select_df['counts_breast'] + select_df['counts_breast'].shift(1)
        select_df['snps_nb_double'] = select_df['counts_nonbreast'] + select_df['counts_nonbreast'].shift(1)
        select_df['start_new'] = select_df['start_region'].shift(1) #.astype(str) #+ "_" + select_df['stop_region'].astype(str)
    #     select_df['stop_new'] = select_df['stop_region']       
        select_df = select_df.reset_index()
        select_df.drop(columns=['index'], inplace=True)
        select_df = select_df.reset_index()
        select_df = select_df[(select_df['index'] % 2 != 0) | (select_df['index'] == select_df['index'].max())]
        select_df['start_new'] = list(select_df['start_new'][:-1])+[list(select_df['start_region'])[-1]]
        select_df.drop(columns=['index'], inplace=True)
    #     select_df['stop_new'] = list(select_df['stop_new'][:-1])+[list(select_df['stop_region'])[-1]]
        df = pd.concat([df, select_df])
    df['start_region'] = df['start_new']
    df.drop(columns=['counts_breast', 'counts_nonbreast', 'start_new'], inplace=True)
    df.rename(columns = {'snps_b_double':'counts_breast', 'snps_nb_double':'counts_nonbreast'}, inplace = True)
    return df


def multiprocess(df, all_snps_b, all_snps_nb, type_df, type_analyse, path_file, select_chrom):
    parts_df_b_nb = np.array_split(df, 20)
    cpus = mp.cpu_count()

    arg_multi_list = []
    for i, df_part in enumerate(parts_df_b_nb):
        arg_multi_list.append((df_part, all_snps_b, all_snps_nb, type_df, type_analyse, path_file, select_chrom))

    pool = Pool(processes=cpus)
    pool.starmap(func=tests.all_test, iterable=arg_multi_list)
    pool.close()
    pool.join()


def all_data(filter_par, path_file, path_db, path_R_b, path_R_nb, select_chrom):
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par, path_file, path_db)
    merge_df = make_merge_df(path_R_b, path_R_nb, select_chrom)
    # df_1000_tests = tests.all_test(merge_df, all_snps_b, all_snps_nb, 'ALL', 'Region_1000', path_file, select_chrom)
    multiprocess(merge_df, all_snps_b, all_snps_nb, 'ALL', 'Region_1000', path_file, select_chrom)
    df_2000 = make_df_2000(merge_df)
    # df_2000_tests = tests.all_test(df_2000, all_snps_b, all_snps_nb, 'ALL', 'Region_2000', path_file, select_chrom)
    multiprocess(df_2000, all_snps_b, all_snps_nb, 'ALL', 'Region_2000', path_file, select_chrom)



def noncoding_data(filter_par, path_file, path_db, path_R_b, path_R_nb, select_chrom):
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(filter_par, path_file, path_db)
    nc_merge_df = make_merge_df(path_R_b, path_R_nb, select_chrom)
    # nc_df_1000_tests = tests.all_test(nc_merge_df, all_snps_b, all_snps_nb, 'NonCoding', 'Region_1000', path_file, select_chrom)
    multiprocess(nc_merge_df, all_snps_b, all_snps_nb, 'NonCoding', 'Region_1000', path_file, select_chrom)
    nc_df_2000 = make_df_2000(nc_merge_df)
    # nc_df_2000_tests = tests.all_test(nc_df_2000, all_snps_b, all_snps_nb, 'NonCoding', 'Region_2000', path_file, select_chrom)
    multiprocess(nc_df_2000, all_snps_b, all_snps_nb, 'NonCoding', 'Region_2000', path_file, select_chrom)


def coding_data(filter_par, path_file, path_db, path_R_b, path_R_nb, select_chrom):
    coding_breast, coding_nonbreast, coding_num_donor_b, coding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_coding_data(filter_par, path_file, path_db)
    c_merge_df = make_merge_df(path_R_b, path_R_nb, select_chrom)
    # c_df_1000_tests = tests.all_test(c_merge_df, all_snps_b, all_snps_nb, 'Coding', 'Region_1000', path_file, select_chrom)
    multiprocess(c_merge_df, all_snps_b, all_snps_nb, 'Coding', 'Region_1000', path_file, select_chrom)
    c_df_2000 = make_df_2000(c_merge_df)
    # c_df_2000_tests = tests.all_test(c_df_2000, all_snps_b, all_snps_nb, 'Coding', 'Region_2000', path_file, select_chrom)
    multiprocess(c_df_2000, all_snps_b, all_snps_nb, 'Coding', 'Region_2000', path_file, select_chrom)



def main():
    config = get_config('calculon')
    path_R = config['path_R']
    path_db = ''
    path_file = config['analyse'] 
    filter_par = False
    select_chrom = sys.argv[1]
    all_data(filter_par, path_file, path_db, f"{path_R}all_breast_ALL.tsv", f"{path_R}all_nonbreast_ALL.tsv", select_chrom)
    noncoding_data(filter_par, path_file, path_db, f"{path_R}noncoding_breast_ALL.tsv", f"{path_R}noncoding_nonbreast_ALL.tsv", select_chrom)
    coding_data(filter_par, path_file, path_db, f"{path_R}coding_breast_ALL.tsv", f"{path_R}coding_nonbreast_ALL.tsv", select_chrom)
    
if __name__ == '__main__':
    main()

