#!/usr/bin/env python3

#Imports
import sys
import numpy as np
import statsmodels.api as sm

sys.path.append(
    '/groups/umcg-wijmenga/tmp04/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')

from config import get_config
from multiprocessing import Pool
import multiprocessing as mp

import get_data as get_data


def make_GT_df(df, type_c, num_donors, select_chrom):
    """
    Creates GT file so that it can be worked well
    :param df: The dataframe with the breast data or the nonbreast data
    :param type_c:  nb: nonbreast or b:breast
    :param num_donors:   Number of donors
    :param select_chrom: The selected chromosome    
    :return:    GT_df_0: File so that we can use the data in follow-up analyses
    """
    select_df = df[['GT2', 'snp_ID', 'chr', 'pos_start', 'pos_end']]
    # Select GT == 1 and GT == 2
    select_1_df = select_df[select_df['GT2'] == 1]
    select_2_df = select_df[select_df['GT2'] == 2]
    # Groupby
    GT1_df = select_1_df.groupby(by=['snp_ID', 'chr', 'pos_start', 'pos_end']).sum().rename(
        columns={'GT2': f'GT_1_{type_c}'})
    GT2_df = select_2_df.groupby(by=['snp_ID', 'chr', 'pos_start', 'pos_end']).sum().rename(
        columns={'GT2': f'GT_2_{type_c}'})
    # Merge GT2_df and GT1_df
    GT_df = GT1_df.merge(GT2_df, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    GT_df_0 = GT_df.fillna(0)
    GT_df_0[f'GT_2_{type_c}'] = GT_df_0[f'GT_2_{type_c}'] / 2
    GT_df_0.sort_values(by=['chr', 'pos_start'], inplace=True)
    GT_df_0[f'GT_0_{type_c}'] = num_donors - (GT_df_0[f'GT_1_{type_c}'] + GT_df_0[f'GT_2_{type_c}'])
    GT_df_0 = GT_df_0[GT_df_0['chr'] == select_chrom]
    return GT_df_0


def cochran_armitage(both_GT, path_file, type_df, select_chrom, i):
    """
    The Cochran–Armitage test for trend is being performed
    :param both_GT: The dataframe with the GT both of breast and nonbreast donors
    :param path_file:  The path to which the files are saved
    :param type_df: Which data is used non-coding, coding or all data
    :param select_chrom:  The selected chromosome
    :param i:        Index of where it is with multiprocess
    :return:    
    """
    both_GT.reset_index(inplace=True)
    p_value_cochran_armitage = list()
    # Loop over rows
    for index, row in both_GT.iterrows():
        # Make continguency table
        contingency_table = np.array(
            ([row['GT_0_b'], row['GT_1_b'], row['GT_2_b']], [row['GT_0_nb'], row['GT_1_nb'], row['GT_2_nb']]))
        table = sm.stats.Table(contingency_table)
        # Perform test
        p_value = table.test_ordinal_association().pvalue
        # Append p-value to list
        p_value_cochran_armitage.append(p_value)
    # Make new column with p-values
    both_GT['p_value_cochran_armitage'] = p_value_cochran_armitage
    both_GT.to_csv(f"{path_file}GT_{type_df}_{select_chrom}_{i}_cochran_armitage.tsv", sep='\t', encoding='utf-8',
                   index=False)


def multiprocess(df, path_file, group, select_chrom):
    """
    Causes function cochran_armitage to be multiprocessed
    :param df: The dataframe with the GT both of breast and nonbreast donors
    :param path_file:  The path to which the files are saved
    :param group:      Which data is used non-coding, coding or all data
    :param select_chrom:     The selected chromosome     
    :return:    
    """
    parts_df_b_nb = np.array_split(df, 20)
    cpus = mp.cpu_count()

    arg_multi_list = []
    for i, df_part in enumerate(parts_df_b_nb):
        arg_multi_list.append((df_part, path_file, group, select_chrom, i))

    pool = Pool(processes=cpus)
    pool.starmap(func=cochran_armitage, iterable=arg_multi_list)
    pool.close()
    pool.join()


def prep_file(breast, num_donor_b, nonbreast, num_donor_nb, non_coding, path_file, select_chrom):
    """
    Prepares the data so that it can be used for the cochran_armitage test
    :param breast:       breast cancer data 
    :param num_donor_b:      The number of breast cancer donors  
    :param nonbreast:      nonbreast cancer data  
    :param num_donor_nb:      The number of nonbreast cancer donors  
    :param non_coding: Which data is used non-coding, coding or all data
    :param path_file:  The path to which the files are saved
    :param select_chrom: The selected chromosome  
    :return:
    """
    # Prepaired the data
    breast_GT = make_GT_df(breast, 'b', num_donor_b, select_chrom)
    nonbreast_GT = make_GT_df(nonbreast, 'nb', num_donor_nb, select_chrom)
    # Merge breast_GT and nonbreast_GT
    both_GT = breast_GT.merge(nonbreast_GT, on=['snp_ID', 'chr', 'pos_start', 'pos_end'], how='outer')
    both_GT['GT_1_b'].fillna(0, inplace=True)
    both_GT['GT_2_b'].fillna(0, inplace=True)
    both_GT['GT_0_b'].fillna(num_donor_b, inplace=True)
    both_GT['GT_1_nb'].fillna(0, inplace=True)
    both_GT['GT_2_nb'].fillna(0, inplace=True)
    both_GT['GT_0_nb'].fillna(num_donor_nb, inplace=True)
    # Call multiprocess
    multiprocess(both_GT, path_file, non_coding, select_chrom)


def main():
    # Call get_config
    config = get_config('calculon')
    # When a database path is listed here, files still need to be created.
    # When this is read, it is assumed that no files need to be created.
    path_db = ''
    path_file = config['analyse']
    filter_par = False
    select_chrom = sys.argv[1].replace('chr', '')
    all_breast, all_nonbreast, all_num_donor_b, all_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(
        filter_par, path_file, path_db)
    prep_file(all_breast, all_num_donor_b, all_nonbreast, all_num_donor_nb, 'ALL', path_file, select_chrom)
    noncoding_breast, noncoding_nonbreast, noncoding_num_donor_b, noncoding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(
        filter_par, path_file, path_db)
    prep_file(noncoding_breast, noncoding_num_donor_b, noncoding_nonbreast, noncoding_num_donor_nb, 'NonCoding',
              path_file, select_chrom)
    coding_breast, coding_nonbreast, coding_num_donor_b, coding_num_donor_nb, all_snps_b, all_snps_nb = get_data.get_coding_data(
        filter_par, path_file, path_db)
    prep_file(coding_breast, coding_num_donor_b, coding_nonbreast, coding_num_donor_nb, 'Coding', path_file,
              select_chrom)


if __name__ == '__main__':
    main()
