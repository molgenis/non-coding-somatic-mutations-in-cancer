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

def make_snp_df(df, type_df, type_c, path_file):
    print('----------')
    df['snp'] = df['chr'].map(str) + '_' + df['pos_start'].map(str) + '_' + df['pos_end'].map(str)
    print(len(list(set(df['donor_ID']))))
    num_donor = len(list(set(df['donor_ID'])))

    df_dict = dict(Counter(list(df['snp'])))
    
    snp_count_df = pd.DataFrame([df_dict.keys(), df_dict.values()]).T
    snp_count_df.columns= ['snp', 'counts']
    max_df = snp_count_df['counts'].max()
    print(snp_count_df['counts'].max())

    snp_count_df[['chr', 'pos_start', 'pos_end']] = snp_count_df['snp'].str.split('_', expand=True)

    snp_count_df.drop('snp', axis=1, inplace=True)

    sort_snp_count_df = snp_count_df.sort_values('counts', ascending=False)

    sort_snp_count_df.rename(columns={'counts': f'counts_{type_c}'}, inplace=True)
    
    df_R = sort_snp_count_df.reset_index().drop([f'counts_{type_c}', 'index'], 1)
    df_R['chr'] ='chr' + df_R['chr'].astype(str)
    df_R.to_csv(f"{path_file}R/{type_df}_{type_c}.tsv", sep='\t', encoding='utf-8', header=None)
    return num_donor, sort_snp_count_df


def run_snp_tests(df_breast, df_nonbreast, type_df, type_analyse, path_file):
    print('\nset snps')
    num_donor_b, sort_snp_count_breast = make_snp_df(df_breast, type_df, 'breast', path_file)
    num_donor_nb, sort_snp_count_nonbreast = make_snp_df(df_nonbreast, type_df, 'nonbreast', path_file)
    
    print('\nmerge dfs')
    sort_snp_count_both = sort_snp_count_breast.merge(sort_snp_count_nonbreast, on=['chr', 'pos_start', 'pos_end'], how='outer')
    sort_snp_count_both_0 = sort_snp_count_both.fillna(0)
    sort_snp_count_both_0.to_csv(f"{path_file}{type_analyse}_{type_df}_both_0.tsv", sep='\t', encoding='utf-8', index=False)
    
    sort_snp_count_both_0 = tests.all_test(sort_snp_count_both_0, num_donor_b, num_donor_nb, type_df, type_analyse, path_file)
    return sort_snp_count_both_0, num_donor_b, num_donor_nb


def all_data(filter_par, path_file, path_db):    
    all_breast, all_nonbreast, num_donor_b, num_donor_nb = get_data.get_all_data(filter_par, path_file, path_db)

    sort_snp_count_all_both_0, all_num_donor_b, all_num_donor_nb = run_snp_tests(all_breast, all_nonbreast, 'ALL', 'per_snp', path_file)
    return sort_snp_count_all_both_0, all_num_donor_b, all_num_donor_nb



def noncoding_data(filter_par, path_file, path_db):
    noncoding_breast, noncoding_nonbreast, num_donor_b, num_donor_nb = get_data.get_noncoding_data(filter_par, path_file, path_db)

    sort_snp_count_noncoding_both_0, noncoding_num_donor_b, noncoding_num_donor_nb = run_snp_tests(noncoding_breast, noncoding_nonbreast, 'NonCoding', 'per_snp', path_file)
    return sort_snp_count_noncoding_both_0, noncoding_num_donor_b, noncoding_num_donor_nb


def coding_data(filter_par, path_file, path_db):
    coding_breast, coding_nonbreast, num_donor_b, num_donor_nb = get_data.get_coding_data(filter_par, path_file, path_db)

    sort_snp_count_coding_both_0, coding_num_donor_b, coding_num_donor_nb = run_snp_tests(coding_breast, coding_nonbreast, 'Coding', 'per_snp', path_file)
    return sort_snp_count_coding_both_0, coding_num_donor_b, coding_num_donor_nb



def main():
    config = get_config()
    path_db = '' #config['database']
    path_file = config['analyse']
    filter_par = False
    sort_snp_count_all_both_0, all_num_donor_b, all_num_donor_nb = all_data(filter_par, path_file, path_db)
    sort_snp_count_noncoding_both_0, noncoding_num_donor_b, noncoding_num_donor_nb = noncoding_data(filter_par, path_file, path_db)
    sort_snp_count_coding_both_0, coding_num_donor_b, coding_num_donor_nb = coding_data(filter_par, path_file, path_db)

if __name__ == '__main__':
    main()