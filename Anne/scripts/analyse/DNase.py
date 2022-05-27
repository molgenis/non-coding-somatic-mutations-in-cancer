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
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

import ast
import get_data as get_data
import tests as tests


def prep_file(path_file, path_save, type_data):
    df = pd.read_csv(path_file, sep='\t')
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
    df_b_nb.to_csv(f"{path_save}b_nb_{type_data}.tsv", sep='\t', encoding='utf-8', index=False)
    return df_b_nb

def run_all(type_data, path_db, path_save):
    filter_par = False
   
    path_file = f"{path_save}{type_data}_chrALL_num_snps.tsv"

    df_b_nb = prep_file(path_file, path_save, type_data)
    print(df_b_nb)
    print(df_b_nb.head())    
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb = get_data.get_all_data(filter_par, path_save, path_db)
    tests_df = tests.all_test(df_b_nb, all_num_donor_b, all_num_donor_nb, 'NonCoding', f'{type_data}', path_save)

def main():
    config = get_config()
    path_db = '' #'D:/Hanze_Groningen/STAGE/lastdb/db_laatste_copy.db' #config['database']
    path_save = config['analyse'] #'D:/Hanze_Groningen/STAGE/UMAP/'
    type_data = 'DNase'
    run_all(type_data, path_db, path_save)



if __name__ == '__main__':
    main()

