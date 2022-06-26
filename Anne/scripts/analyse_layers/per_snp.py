#!/usr/bin/env python3

#Imports
import sys
import pandas as pd
from collections import Counter

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
import get_data as get_data
import tests as tests


def make_snp_df(df, type_df, type_c, path_file):
    """
    Creates the snp dataframes to work with
    :param df: The dataframe with the breast cancer data or the nonbreast cancer data
    :param type_df:  type of dataframe (all data, noncoding data or coding data)
    :param type_c:  breast cancer or nonbreast cancer
    :param path_file:  Path after which files are saved
    :return:    num_donor: Number of donors
                sort_snp_count_df: Data frame with 1 column with SNP IDs and the other column with how
                                   often this snp occurs
    """

    df['snp'] = df['chr'].map(str) + '_' + df['pos_start'].map(str) + '_' + df['pos_end'].map(str)
    # Number of donors
    num_donor = len(list(set(df['donor_ID'])))
    # Dictionary with the key snp ID and the value how often this ID occurs
    df_dict = dict(Counter(list(df['snp'])))
    # Make dataframe of dictionary
    snp_count_df = pd.DataFrame([df_dict.keys(), df_dict.values()]).T
    snp_count_df.columns = ['snp', 'counts']
    snp_count_df[['chr', 'pos_start', 'pos_end']] = snp_count_df['snp'].str.split('_', expand=True)
    snp_count_df.drop('snp', axis=1, inplace=True)
    # Sort dataframe on counts
    sort_snp_count_df = snp_count_df.sort_values('counts', ascending=False)
    sort_snp_count_df.rename(columns={'counts': f'counts_{type_c}'}, inplace=True)
    # For in R
    df_R = sort_snp_count_df.reset_index().drop([f'counts_{type_c}', 'index'], 1)
    df_R['chr'] = 'chr' + df_R['chr'].astype(str)
    df_R.to_csv(f"{path_file}R/{type_df}_{type_c}.tsv", sep='\t', encoding='utf-8', header=None)
    return num_donor, sort_snp_count_df


def run_snp_tests(df_breast, df_nonbreast, type_df, type_analyse, path_file, select_chrom, i):
    """
    Prepare the files and ensure that all tests are run multiprocessed
    :param df_breast: Dataframe with breast cancer data
    :param df_nonbreast:  Dataframe with nonbreast cancer data
    :param type_df: type of dataframe (all data, noncoding data or coding data)
    :param type_analyse: What kind of analysis is being performed. (in this case it's per_snp every time)
    :param path_file:  Path after which files are saved
    :param select_chrom: The selected chromosome  (chr0)
    :param i: always 0   
    :return:    
    """
    # Make dataframes for breast cancer and nonbreast cancer
    num_donor_b, sort_snp_count_breast = make_snp_df(df_breast, type_df, 'breast', path_file)
    num_donor_nb, sort_snp_count_nonbreast = make_snp_df(df_nonbreast, type_df, 'nonbreast', path_file)
    # Merge sort_snp_count_breast and sort_snp_count_nonbreast
    sort_snp_count_both = sort_snp_count_breast.merge(sort_snp_count_nonbreast, on=['chr', 'pos_start', 'pos_end'],
                                                      how='outer')
    sort_snp_count_both_0 = sort_snp_count_both.fillna(0)
    sort_snp_count_both_0.to_csv(f"{path_file}{type_analyse}_{type_df}_both_0.tsv", sep='\t', encoding='utf-8',
                                 index=False)
    # Run the tests
    sort_snp_count_both_0 = tests.all_test(sort_snp_count_both_0, num_donor_b, num_donor_nb, type_df, type_analyse,
                                           path_file, select_chrom, i)


def main():
    # Call get_config
    config = get_config('gearshift')
    # When a database path is listed here, files still need to be created.
    # When this is read, it is assumed that no files need to be created.
    path_db = ''
    path_file = config['analyse']
    filter_par = False
    select_chrom = 'chr0'
    i = '0'
    # Run test over all data, noncoding data and coding data
    all_breast, all_nonbreast, num_donor_b, num_donor_nb, all_snps_b, all_snps_nb = get_data.get_all_data(filter_par,
                                                                                                          path_file,
                                                                                                          path_db)
    run_snp_tests(all_breast, all_nonbreast, 'ALL', 'per_snp', path_file, select_chrom, i)
    noncoding_breast, noncoding_nonbreast, num_donor_b, num_donor_nb, all_snps_b, all_snps_nb = get_data.get_noncoding_data(
        filter_par, path_file, path_db)
    run_snp_tests(noncoding_breast, noncoding_nonbreast, 'NonCoding', 'per_snp', path_file, select_chrom, i)
    coding_breast, coding_nonbreast, num_donor_b, num_donor_nb, all_snps_b, all_snps_nb = get_data.get_coding_data(
        filter_par, path_file, path_db)
    run_snp_tests(coding_breast, coding_nonbreast, 'Coding', 'per_snp', path_file, select_chrom, i)


if __name__ == '__main__':
    main()
