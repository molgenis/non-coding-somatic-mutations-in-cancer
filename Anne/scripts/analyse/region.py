import sys
# import multiprocessing as mp
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats.distributions import chi2
from bioinfokit import analys, visuz
from scipy.stats import fisher_exact
import time
from scipy.special import factorial
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from fisher import pvalue_npy
from scipy.stats import chi2_contingency
from scipy.stats import uniform, randint
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

import get_data as get_data
import tests as tests
import per_snp as per_snp

def make_merge_df(path_breast, path_nonbreast):
    colnames=['index', 'counts_breast', 'chr', 'start_region', 'stop_region']
    breast = pd.read_csv(path_breast, sep='\t', header=None, names=colnames)
    breast.sort_values(['chr', 'start_region'], inplace=True)
    
    colnames=['index', 'counts_nonbreast', 'chr', 'start_region', 'stop_region']
    nonbreast = pd.read_csv(path_nonbreast, sep='\t', header=None, names=colnames)   
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



def all_data(filter_par, path_file, path_db, path_R_b, path_R_nb):
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb = get_data.get_all_data(filter_par, path_file, path_db)
    merge_df = make_merge_df(path_R_b, path_R_nb)
    df_1000_tests = tests.all_test(merge_df, all_num_donor_b, all_num_donor_nb, 'ALL', 'Region_1000', path_file)
    df_2000 = make_df_2000(merge_df)
    df_2000_tests = tests.all_test(df_2000, all_num_donor_b, all_num_donor_nb, 'ALL', 'Region_2000', path_file)


def noncoding_data(filter_par, path_file, path_db, path_R_b, path_R_nb):
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb = get_data.get_noncoding_data(filter_par, path_file, path_db)
    nc_merge_df = make_merge_df(path_R_b, path_R_nb)
    nc_df_1000_tests = tests.all_test(nc_merge_df, noncoding_num_donor_b, noncoding_num_donor_nb, 'NonCoding', 'Region_1000', path_file)
    nc_df_2000 = make_df_2000(nc_merge_df)
    nc_df_2000_tests = tests.all_test(nc_df_2000, noncoding_num_donor_b, noncoding_num_donor_nb, 'NonCoding', 'Region_2000', path_file)

def coding_data(filter_par, path_file, path_db, path_R_b, path_R_nb):
    coding_breast, coding_nonbreast, coding_num_donor_b, coding_num_donor_nb = get_data.get_coding_data(filter_par, path_file, path_db)
    c_merge_df = make_merge_df(path_R_b, path_R_nb)
    c_df_1000_tests = tests.all_test(c_merge_df, coding_num_donor_b, coding_num_donor_nb, 'Coding', 'Region_1000', path_file)
    c_df_2000 = make_df_2000(c_merge_df)
    c_df_2000_tests = tests.all_test(c_df_2000, coding_num_donor_b, coding_num_donor_nb, 'Coding', 'Region_2000', path_file)


def main():
    path_R = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/analyse/R/' #'D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/all/'
    config = get_config()
    path_db = '' #'D:/Hanze_Groningen/STAGE/lastdb/db_laatste_copy.db' #config['database']
    path_file = config['analyse'] #config['analyse'] #'D:/Hanze_Groningen/STAGE/lastdb/'
    filter_par = False
    all_data(filter_par, path_file, path_db, f"{path_R}all_breast_ALL.tsv", f"{path_R}all_nonbreast_ALL.tsv")
    noncoding_data(filter_par, path_file, path_db, f"{path_R}noncoding_breast_ALL.tsv", f"{path_R}noncoding_nonbreast_ALL.tsv")
    coding_data(filter_par, path_file, path_db, f"{path_R}coding_breast_ALL.tsv", f"{path_R}coding_nonbreast_ALL.tsv")
    
if __name__ == '__main__':
    main()

