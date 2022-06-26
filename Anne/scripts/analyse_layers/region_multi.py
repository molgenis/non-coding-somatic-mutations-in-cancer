#!/usr/bin/env python3

#Imports
import sys
import pandas as pd
import numpy as np

sys.path.append(
    '/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from multiprocessing import Pool
import multiprocessing as mp
import get_data as get_data
import tests as tests


def make_merge_df(path_breast, path_nonbreast, select_chrom):
    """
    Merged the two data frames to create a data frame with both breast cancer counts and non breast cancer counts
    :param path_breast: Path to R files breast
    :param path_nonbreast:  Path to R files nonbreast
    :param select_chrom:       The selected chromosome   
    :return:   merge_df:       1000 bp data frame with breast cancer and non breast cancer counts 
    """
    colnames = ['index', 'counts_breast', 'chr', 'start_region', 'stop_region']
    breast = pd.read_csv(path_breast, sep='\t', header=None, names=colnames)
    breast = breast[breast['chr'] == select_chrom]
    breast.sort_values(['chr', 'start_region'], inplace=True)

    colnames = ['index', 'counts_nonbreast', 'chr', 'start_region', 'stop_region']
    nonbreast = pd.read_csv(path_nonbreast, sep='\t', header=None, names=colnames)
    nonbreast = nonbreast[nonbreast['chr'] == select_chrom]
    nonbreast.sort_values(['chr', 'start_region'], inplace=True)

    merged_df = breast.merge(nonbreast, on=['chr', 'start_region', 'stop_region'], how='outer')
    merged_df.drop(columns=['index_x', 'index_y'], inplace=True)
    return merged_df


def make_df_2000(merge_df):
    """
    Creates the 2000 bp data frame from the 1000 bp data frame
    :param merge_df:       1000 bp data frame with breast cancer and non breast cancer counts
    :return: df:   2000 bp data frame with breast cancer and non breast cancer counts
    """
    df = pd.DataFrame(
        columns=['counts_breast', 'chr', 'start_region', 'stop_region', 'counts_nonbreast', 'snps_b_double',
                 'snps_nb_double'])
    # Loop over chromosomes
    for i in list(set(merge_df['chr'])):
        # Filter dataframe on chromosome
        select_df = merge_df[merge_df['chr'] == i]
        # Slide all data
        select_df['snps_b_double'] = select_df['counts_breast'] + select_df['counts_breast'].shift(1)
        select_df['snps_nb_double'] = select_df['counts_nonbreast'] + select_df['counts_nonbreast'].shift(1)
        select_df['start_new'] = select_df['start_region'].shift(1)
        select_df = select_df.reset_index()
        select_df.drop(columns=['index'], inplace=True)
        select_df = select_df.reset_index()
        select_df = select_df[(select_df['index'] % 2 != 0) | (select_df['index'] == select_df['index'].max())]
        select_df['start_new'] = list(select_df['start_new'][:-1]) + [list(select_df['start_region'])[-1]]
        select_df.drop(columns=['index'], inplace=True)
        df = pd.concat([df, select_df])
    df['start_region'] = df['start_new']
    df.drop(columns=['counts_breast', 'counts_nonbreast', 'start_new'], inplace=True)
    df.rename(columns={'snps_b_double': 'counts_breast', 'snps_nb_double': 'counts_nonbreast'}, inplace=True)
    return df


def multiprocess(df, all_snps_b, all_snps_nb, type_df, type_analyse, path_file, select_chrom):
    """
    Multiprocess the tests
    :param df: The dataframe with the info
    :param all_snps_b:  SNPs for breast cancer
    :param all_snps_nb:  SNPs for nonbreast cancer
    :param type_df:  noncoding, coding or mix 
    :param type_analyse:  What analysis is done in this case is either 1000 bp regions or 2000 bp regions
    :param path_file:  Path after which files are saved
    :param select_chrom: The selected chromosome     
    :return:    
    """
    parts_df_b_nb = np.array_split(df, 20)
    cpus = mp.cpu_count()

    arg_multi_list = []
    for i, df_part in enumerate(parts_df_b_nb):
        arg_multi_list.append((df_part, all_snps_b, all_snps_nb, type_df, type_analyse, path_file, select_chrom))

    pool = Pool(processes=cpus)
    pool.starmap(func=tests.all_test, iterable=arg_multi_list)
    pool.close()
    pool.join()


def get_all_data(snps_b, snps_nb, path_file, path_R_b, path_R_nb, select_chrom, type_data):
    """
    Make sure the right dataframe is formed and run multiprocess
    :param snps_b: SNPs for breast cancer
    :param snps_nb: SNPs for nonbreast cancer
    :param path_file:  Path after which files are saved
    :param path_db: Path to the database
    :param path_R_b: Path to R files breast
    :param path_R_nb:  Path to R files nonbreast
    :param select_chrom:      The selected chromosome    
    :param type_data: noncoding, coding or mix   
    :return:    
    """
    # Run for 1000 bp
    merge_df = make_merge_df(path_R_b, path_R_nb, select_chrom)
    multiprocess(merge_df, snps_b, snps_nb, type_data, 'Region_1000', path_file, select_chrom)
    # Run for 2000 bp
    df_2000 = make_df_2000(merge_df)
    multiprocess(df_2000, snps_b, snps_nb, type_data, 'Region_2000', path_file, select_chrom)


def main():
    # Call get_config
    config = get_config('calculon')
    path_R = config['path_R']
    # When a database path is listed here, files still need to be created.
    # When this is read, it is assumed that no files need to be created.
    path_db = ''
    path_file = config['analyse']
    filter_par = False
    select_chrom = sys.argv[1]
    # Run for all data, noncoding data and coding data
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(
        filter_par, path_file, path_db)
    get_all_data(all_snps_b, all_snps_nb, path_file, f"{path_R}all_breast_ALL.tsv", f"{path_R}all_nonbreast_ALL.tsv",
                 select_chrom, 'ALL')
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(
        filter_par, path_file, path_db)
    get_all_data(all_snps_b, all_snps_nb, path_file, f"{path_R}noncoding_breast_ALL.tsv",
                 f"{path_R}noncoding_nonbreast_ALL.tsv", select_chrom, 'NonCoding')
    coding_breast, coding_nonbreast, coding_num_donor_b, coding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_coding_data(
        filter_par, path_file, path_db)
    get_all_data(all_snps_b, all_snps_nb, path_file, f"{path_R}coding_breast_ALL.tsv",
                 f"{path_R}coding_nonbreast_ALL.tsv", select_chrom, 'Coding')


if __name__ == '__main__':
    main()
