from Database import Database
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
# sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
# from config import get_config


def get_data_db(filter_par, path_file, path_db):
    """
    
    """
    if filter_par:
        # Database connection
        db = Database(path_db)
        df_whole = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                                    sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, sum_dosage_GT.GT2,
                                    snp.chr, snp.pos_start, snp.pos_end, snp.ref, snp.alt,
                                    snp.UCNE, snp.TFBS, snp.DNase, 
                                    snp.in_transcript, snp.in_coding, snp.in_exon,
                                    snp.before_gene, snp.after_gene, 
                                    snp.ID_eQTL, snp.eQTL, snp.close_eQTL_3000
                            FROM project, sum_dosage_GT, snp 
                            WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                                    sum_dosage_GT.donor_project_ID = project.ID AND 
                                    (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                                    sum_dosage_GT.total_read_count_sum >= 33;''', db.mydb_connection)
        db.close()
        df_whole.to_csv(f"{path_file}df_whole.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        df_whole = pd.read_csv(f"{path_file}df_whole.tsv", sep='\t')

    df_whole = df_whole.loc[df_whole['chr'] != 'MT']
    return df_whole

def filter_whole_cancer(filter_par, df_whole, path_file):
    """
    
    """
    if filter_par:
        all_breast =  df_whole.loc[df_whole['cancer'] == 'Breast']
        all_breast.to_csv(f"{path_file}all_breast.tsv", sep='\t', encoding='utf-8', index=False)
        all_nonbreast = df_whole.loc[df_whole['cancer'] != 'Breast']
        all_nonbreast.to_csv(f"{path_file}all_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        all_breast = pd.read_csv(f"{path_file}all_breast.tsv", sep='\t')
        all_nonbreast = pd.read_csv(f"{path_file}all_nonbreast.tsv", sep='\t')
    return all_breast, all_nonbreast


def filter_noncoding(filter_par, df_whole, path_file):
    """
    
    """
    if filter_par:
        noncoding_df = df_whole.loc[(df_whole['in_transcript'] == 0) & (df_whole['in_coding'] == 0) & (df_whole['in_exon'] == 0)]
        noncoding_df.to_csv(f"{path_file}noncoding_df.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        noncoding_df = pd.read_csv(f"{path_file}noncoding_df.tsv", sep='\t')
    return noncoding_df


def filter_noncoding_cancer(filter_par, noncoding_df, path_file):
    """
    
    """
    if filter_par:
        noncoding_breast = noncoding_df.loc[noncoding_df['cancer'] == 'Breast']
        noncoding_breast.to_csv(f"{path_file}noncoding_breast.tsv", sep='\t', encoding='utf-8', index=False)
        noncoding_nonbreast = noncoding_df.loc[noncoding_df['cancer'] != 'Breast']
        noncoding_nonbreast.to_csv(f"{path_file}noncoding_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        noncoding_breast = pd.read_csv(f"{path_file}noncoding_breast.tsv", sep='\t')
        noncoding_nonbreast = pd.read_csv(f"{path_file}noncoding_nonbreast.tsv", sep='\t')
    return noncoding_breast, noncoding_nonbreast




def filter_coding(filter_par, df_whole, path_file):
    """
    
    """
    if filter_par:
        coding_df = df_whole.loc[(df_whole['in_transcript'] == 1) | (df_whole['in_coding'] == 1) | (df_whole['in_exon'] == 1)]
        coding_df.to_csv(f"{path_file}coding_df.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        coding_df = pd.read_csv(f"{path_file}coding_df.tsv", sep='\t')
    return coding_df

def filter_coding_cancer(filter_par, coding_df, path_file):
    """
    
    """
    if filter_par:
        coding_breast = coding_df.loc[coding_df['cancer'] == 'Breast']
        coding_breast.to_csv(f"{path_file}coding_breast.tsv", sep='\t', encoding='utf-8', index=False)
        coding_nonbreast = coding_df.loc[coding_df['cancer'] != 'Breast']
        coding_nonbreast.to_csv(f"{path_file}coding_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        coding_breast = pd.read_csv(f"{path_file}coding_breast.tsv", sep='\t')
        coding_nonbreast = pd.read_csv(f"{path_file}coding_nonbreast.tsv", sep='\t')
    return coding_breast, coding_nonbreast


def get_all_data(filter_par, path_file, path_db):
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        all_breast, all_nonbreast = filter_whole_cancer(filter_par, df_whole, path_file)
    else:
        all_breast, all_nonbreast = filter_whole_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(all_breast['donor_ID'])))
    num_donor_nb = len(list(set(all_nonbreast['donor_ID'])))

    return all_breast, all_nonbreast, num_donor_b, num_donor_nb

def get_noncoding_data(filter_par, path_file, path_db):
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        noncoding_df = filter_noncoding(filter_par, df_whole, path_file)
        noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(filter_par, noncoding_df, path_file)
    else:
        noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(noncoding_breast['donor_ID'])))
    num_donor_nb = len(list(set(noncoding_nonbreast['donor_ID'])))

    return noncoding_breast, noncoding_nonbreast, num_donor_b, num_donor_nb

def get_coding_data(filter_par, path_file, path_db):
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        coding_df = filter_coding(filter_par, df_whole, path_file)
        coding_breast, coding_nonbreast = filter_coding_cancer(filter_par, coding_df, path_file)
    else:
        coding_breast, coding_nonbreast = filter_coding_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(coding_breast['donor_ID'])))
    num_donor_nb = len(list(set(coding_nonbreast['donor_ID'])))

    return coding_breast, coding_nonbreast, num_donor_b, num_donor_nb

def main():
    config = get_config()
    path_file = config['analyse']
    path_db = config['database_get_data']
    df_whole = get_data_db(True, path_file, path_db)
    all_breast, all_nonbreast = filter_whole_cancer(True, df_whole, path_file)
    noncoding_df = filter_noncoding(True, df_whole, path_file)
    noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(True, noncoding_df, path_file)
    coding_df = filter_coding(True, df_whole, path_file)
    coding_breast, coding_nonbreast = filter_coding_cancer(True, coding_df, path_file)

if __name__ == '__main__':
    main()
