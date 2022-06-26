#!/usr/bin/env python3

#Imports
from statsmodels.sandbox.stats.multicomp import multipletests
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database


def run_multiple_testing_correction(df, alpha, method, type_test):
    """
    Performs the bonferroni test or benjamini hochberg
    :param df:     The data frame with the p-values
    :param alpha:  The alpha being used. When the p-value is less than alpha it is significant
    :param method: Method for running (bonferroni or bh)
    :param type_test:    Which test manual chi-square, chi-square or fisher exactly test
    :return:    df: df with a new column
    """
    # Selects a column and makes it a list
    list_p_values = list(df[f'p_value_{type_test}'])
    # Performs one of two tests
    rej, p_adjusted = multipletests(list_p_values, alpha=alpha, method=method)[:2]
    # Add new column
    df[f'{method}_{type_test}'] = p_adjusted
    return df


def bonferroni_correction(df, alpha, method, type_test):
    """
    Performs the bonferroni correction
    :param df:     The data frame with the p-values
    :param alpha:  The alpha being used. When the p-value is less than alpha it is significant
    :param method:  Method for running (bonferroni or bh)
    :param type_test:    Which test manual chi-square, chi-square or fisher exactly test 
    :return:    df: df with a new column
    """
    df[f'{method}_{type_test}'] = df[f'p_value_{type_test}'] < (alpha / len(set(df)))
    return df


def benjamini_hochberg_correction(df, alpha, type_test):
    """
    Performs the benjamini hochberg correctie
    :param df:        The data frame with the p-values
    :param alpha:     The alpha being used. When the p-value is less than alpha it is significant
    :param type_test:  Which test manual chi-square, chi-square or fisher exactly test      
    :return:    
    """
    # sort p-values
    df = df.sort_values(f'p_value_{type_test}').reset_index(drop=True)
    # assign a rank
    df[f'rank_{type_test}'] = df.index + 1
    # Compute BH-critical value
    """
    (i/m)*a
    i = the individual p-value’s rank,
    m = total number of tests,
    alpha = the false discovery rate (a percentage, chosen by you).
    """
    df[f'(i/m)a_{type_test}'] = (df[f'rank_{type_test}'] / len(df)) * alpha
    # Find the largest rank for which the p-value is less than the corresponding critical value
    rank_begin = 0
    list_true_false = list()
    for index, row in df.iterrows():
        if (row[f'p_value_{type_test}'] < row[f'(i/m)a_{type_test}']) and row[f'rank_{type_test}'] == (rank_begin + 1):
            rank_begin = row[f'rank_{type_test}']
            list_true_false.append(True)
        else:
            list_true_false.append(False)
    df[f'check_bh_{type_test}'] = list_true_false
    return df


def get_significant(df, type_test, alpha, RR):
    """
    Checks which values are all significant
    :param df:     The data frame with the p-values
    :param alpha:  The alpha being used. When the p-value is less than alpha it is significant
    :param type_test:  Which test manual chi-square, chi-square or fisher exactly test  
    :param RR:      True if it contains the relative risk 
    :return:    genes_sig: Significant regions at a p-value on which no multiple correction has been performed
                bon_sig: significant regions at a p-value with bonferroni correction
                bh_sig: significant genes at a p-value run on benjamini hochberg
    """
    if RR:
        print(f'---------RR--------')
        df_select = df[['info', f'{type_test}']]
        RR_p = df_select[df_select[f'{type_test}'] == True]
        RR_sig = RR_p['info']
        print(f'RR: {len(RR_sig)}')
        return RR_sig
    else:
        df_select = df[['info', f'p_value_{type_test}', f'bonferroni_{type_test}', f'check_bh_{type_test}']]
        print(f'---------{type_test}--------')
        normal_p = df_select[df_select[f'p_value_{type_test}'] <= alpha]
        genes_sig = normal_p['info']
        print(f'nomal: {len(genes_sig)}')

        bon_p = df_select[df_select[f'bonferroni_{type_test}'] == True]
        bon_sig = bon_p['info']
        print(f'bonferroni: {len(bon_sig)}')

        bh_p = df_select[df_select[f'check_bh_{type_test}'] == True]
        bh_sig = bh_p['info']
        print(f'check_bh: {len(bh_sig)}')
        return genes_sig, bon_sig, bh_sig


def get_overlap(self_X2_df, X2_df, F_df):
    """
    Get the overlap between the three different tests in terms of genes
    :param self_X2_df: A list of significant regions from the manual chi-square
    :param X2_df:      A list of significant regions from the chi-square
    :param F_df:       A list of significant regions from the fisher exact test
    :return:    elements_in_all: regions that show up significantly in all tests
    """
    print('---ALL')
    all_values = set(self_X2_df) | set(X2_df) | set(F_df)
    print(f'All genes ({len(all_values)})')

    elements_in_all = list(set.intersection(*map(set, [self_X2_df, X2_df, F_df])))
    print(f'Gene in all list ({len(elements_in_all)})')

    # selfX2 and X2
    print('---selfX2 and X2')
    in_selfX2_X2 = set(self_X2_df) & set(X2_df)
    print(f'in selfX2 and X2 ({len(in_selfX2_X2)})')

    dif_selfX2_X2 = set(self_X2_df) - set(X2_df)
    print(f'in selfX2 NOT in X2 ({len(dif_selfX2_X2)})')

    dif_X2_selfX2 = set(X2_df) - set(self_X2_df)
    print(f'in X2 NOT in selfx2 ({len(dif_X2_selfX2)})')

    # selfX2 and F
    print('---selfX2 and F')
    in_selfX2_F = set(self_X2_df) & set(F_df)
    print(f'in selfX2 and F ({len(in_selfX2_F)})')

    dif_selfX2_F = set(self_X2_df) - set(F_df)
    print(f'in selfX2 NOT in F ({len(dif_selfX2_F)})')

    dif_F_selfX2 = set(F_df) - set(self_X2_df)
    print(f'in F NOT in selfX2 ({len(dif_F_selfX2)})')

    # F and X2
    print('---F and X2')
    in_F_X2 = set(F_df) & set(X2_df)
    print(f'in F and X2 ({len(in_F_X2)})')

    dif_F_X2 = set(F_df) - set(X2_df)
    print(f'in F NOT in X2 ({len(dif_F_X2)})')

    dif_X2_F = set(X2_df) - set(F_df)
    print(f'in X2 NOT in F ({len(dif_X2_F)})\n')
    return elements_in_all


def search(df_select, type_file, non_coding, path_analyse, GT, fc):
    """
    Calls the functions to perform the multiple testing corrections and to determine the significant
    regions and their overlap
    :param df_select:  The dataframe with the p-values
    :param type_file:  What type it is, TFBS, GT or any of the other options
    :param non_coding:  Whether they are non-coding areas or coding or a mix of both
    :param path_analyse:  The path where the new files will be saved
    :param GT:  If it is GT then it is true otherwise false
    :param fc:      Which foldchange is used
    :return: elements_in_all_normal/cochran_armitage_sig_normal: Significant regions occurring with all normal p-values
             elements_in_all_bon/cochran_armitage_bon_sig: Significant regions occurring with all bonferroni p-values
             elements_in_all_bh/cochran_armitage_bh_sig: Significant regions occurring with all benjamini hochberg p-values
             elements_snps_all_MTC/elements_snps_all_MTC: Significant regions occurring in both the normal p-values,
                                                         the bonferroni p-values and the benjamini hochberg p-values
    """
    alpha = 0.05
    if GT:
        # Perform the two multiple testing correction
        method = 'bonferroni'
        df_select = bonferroni_correction(df_select, alpha, method, 'cochran_armitage')
        method = 'fdr_bh'
        df_select = benjamini_hochberg_correction(df_select, alpha, 'cochran_armitage')
        # Save the file
        df_select.to_csv(f"{path_analyse}{type_file}_{non_coding}_{fc}_MTC.tsv", sep='\t', encoding='utf-8',
                         index=False)
        RR = False
        # Get the significant regions
        cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig = get_significant(df_select,
                                                                                                         'cochran_armitage',
                                                                                                         alpha, RR)
        elements_snps_all_MTC = list(set.intersection(
            *map(set, [cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig])))
        return cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig, elements_snps_all_MTC
    else:
        # Perform the two multiple testing correction
        method = 'bonferroni'
        df_select = bonferroni_correction(df_select, alpha, method, 'X2_self')
        df_select = bonferroni_correction(df_select, alpha, method, 'X2')
        df_select = bonferroni_correction(df_select, alpha, method, 'F')

        method = 'fdr_bh'
        df_select = benjamini_hochberg_correction(df_select, alpha, 'X2_self')
        df_select = benjamini_hochberg_correction(df_select, alpha, 'X2')
        df_select = benjamini_hochberg_correction(df_select, alpha, 'F')
        # Get the significant regions
        df_select.to_csv(f"{path_analyse}{type_file}_{non_coding}_{fc}_MTC.tsv", sep='\t', encoding='utf-8',
                         index=False)
        RR = False
        X2_self_sig_normal, X2_self_bon_sig, X2_self_bh_sig = get_significant(df_select, 'X2_self', alpha, RR)
        X2_sig_normal, X2_bon_sig, X2_bh_sig = get_significant(df_select, 'X2', alpha, RR)
        F_sig_normal, F_bon_sig, F_bh_sig = get_significant(df_select, 'F', alpha, RR)
        RR = True
        RR_sig = get_significant(df_select, 'one_in_interval_bon', alpha, RR)

        # Get the overlap between the significant regions
        print('\nNORMAL')
        elements_in_all_normal = get_overlap(X2_self_sig_normal, X2_sig_normal, F_sig_normal)
        print('\nBON')
        elements_in_all_bon = get_overlap(X2_self_bon_sig, X2_bon_sig, F_bon_sig)
        print('\nBH')
        elements_in_all_bh = get_overlap(X2_self_bh_sig, X2_bh_sig, F_bh_sig)
        elements_snps_all_MTC = list(
            set.intersection(*map(set, [elements_in_all_normal, elements_in_all_bon, elements_in_all_bh])))
        return elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC
