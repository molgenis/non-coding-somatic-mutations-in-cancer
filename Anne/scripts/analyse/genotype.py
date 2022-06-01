import sys
# import multiprocessing as mp
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import fisher_exact
import time
from scipy.special import factorial
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency
from scipy.stats import uniform, randint
import statsmodels.api as sm
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

import get_data as get_data

def make_GT_df(df, type_c, num_donors):
    select_df = df[['GT2', 'snp_ID', 'chr', 'pos_start', 'pos_end']]
    select_1_df = select_df[select_df['GT2'] == 1]
    select_2_df = select_df[select_df['GT2'] == 2]
    GT1_df = select_1_df.groupby(by=['snp_ID', 'chr', 'pos_start', 'pos_end']).sum().rename(columns = {'GT2':f'GT_1_{type_c}'})
    GT2_df = select_2_df.groupby(by=['snp_ID', 'chr', 'pos_start', 'pos_end']).sum().rename(columns = {'GT2':f'GT_2_{type_c}'})
    GT_df = GT1_df.merge(GT2_df, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    GT_df_0 = GT_df.fillna(0)
    GT_df_0[f'GT_2_{type_c}'] = GT_df_0[f'GT_2_{type_c}'] /2
    GT_df_0.sort_values(by=['chr', 'pos_start'], inplace = True)
    GT_df_0[f'GT_0_{type_c}'] = num_donors - (GT_df_0[f'GT_1_{type_c}']+GT_df_0[f'GT_2_{type_c}'])
    return GT_df_0

def cochran_armitage(both_GT, path_file, type_df):
    both_GT.reset_index(inplace=True)
    p_value_cochran_armitage = list()
    for index, row in both_GT.iterrows():
        # pd.crosstab([a,d], [b, c, ], rownames=['breast', 'nonbreast'], colnames=['0', '1', '2'])
        contingency_table = np.array(([row['GT_0_b'], row['GT_1_b'], row['GT_2_b']], [row['GT_0_nb'], row['GT_1_nb'], row['GT_2_nb']]))
        print(contingency_table)
        table = sm.stats.Table(contingency_table)
        print(table)
        p_value = table.test_ordinal_association().pvalue #row_scores=np.array([1,0]), col_scores=np.array([0,1,2])
        print(p_value)
        p_value_cochran_armitage.append(p_value)
    both_GT['p_value_cochran_armitage'] = p_value_cochran_armitage
    both_GT.to_csv(f"{path_file}GT_{type_df}_cochran_armitage.tsv", sep='\t', encoding='utf-8', index=False)
    return both_GT


def all_data(filter_par, path_file, path_db):
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par, path_file, path_db)
    breast_GT = make_GT_df(all_breast, 'b', all_num_donor_b)
    nonbreast_GT = make_GT_df(all_nonbreast, 'nb', all_num_donor_nb)
    both_GT = breast_GT.merge(nonbreast_GT, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    both_GT['GT_1_b'].fillna(0,inplace=True)
    both_GT['GT_2_b'].fillna(0,inplace=True)
    both_GT['GT_0_b'].fillna(all_num_donor_b,inplace=True)
    both_GT['GT_1_nb'].fillna(0,inplace=True)
    both_GT['GT_2_nb'].fillna(0,inplace=True)
    both_GT['GT_0_nb'].fillna(all_num_donor_nb,inplace=True)
    both_GT = cochran_armitage(both_GT, path_file, 'ALL')



def noncoding_data(filter_par, path_file, path_db):
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(filter_par, path_file, path_db)
    breast_GT = make_GT_df(noncoding_breast, 'b', noncoding_num_donor_b)
    nonbreast_GT = make_GT_df(noncoding_nonbreast, 'nb', noncoding_num_donor_nb)
    both_GT = breast_GT.merge(nonbreast_GT, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    both_GT['GT_1_b'].fillna(0,inplace=True)
    both_GT['GT_2_b'].fillna(0,inplace=True)
    both_GT['GT_0_b'].fillna(noncoding_num_donor_b,inplace=True)
    both_GT['GT_1_nb'].fillna(0,inplace=True)
    both_GT['GT_2_nb'].fillna(0,inplace=True)
    both_GT['GT_0_nb'].fillna(noncoding_num_donor_nb,inplace=True)       
    both_GT = cochran_armitage(both_GT, path_file, 'NonCoding')
    



def coding_data(filter_par, path_file, path_db):
    coding_breast, coding_nonbreast, coding_num_donor_b, coding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_coding_data(filter_par, path_file, path_db)
    breast_GT = make_GT_df(coding_breast, 'b', coding_num_donor_b)
    nonbreast_GT = make_GT_df(coding_nonbreast, 'nb', coding_num_donor_nb)
    both_GT = breast_GT.merge(nonbreast_GT, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    both_GT['GT_1_b'].fillna(0,inplace=True)
    both_GT['GT_2_b'].fillna(0,inplace=True)
    both_GT['GT_0_b'].fillna(coding_num_donor_b,inplace=True)
    both_GT['GT_1_nb'].fillna(0,inplace=True)
    both_GT['GT_2_nb'].fillna(0,inplace=True)
    both_GT['GT_0_nb'].fillna(coding_num_donor_nb,inplace=True)
    both_GT = cochran_armitage(both_GT, path_file, 'Coding')



def main():
    config = get_config('gearshift')
    path_db = '' 
    path_file = config['analyse'] #config['analyse'] 'D:/Hanze_Groningen/STAGE/lastdb/'
    filter_par = False
    all_data(filter_par, path_file, path_db)
    noncoding_data(filter_par, path_file, path_db)
    coding_data(filter_par, path_file, path_db)


if __name__ == '__main__':
    main()
