import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from Database import Database


def bonferroni_correction(df, alpha, method, type_test):
    df[f'{method}_{type_test}'] = df[f'p_value_{type_test}'] < (alpha/len(set(df)))
    return df

def benjamini_hochberg_correction(df, alpha, type_test):
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
    df[f'(i/m)a_{type_test}'] = (df[f'rank_{type_test}']/len(df))*alpha
    # Find the largest rank for which the p-value is less than the corresponding critical value
    rank_begin = 0
    list_true_false = list()
    for index, row in df.iterrows():
        if (row[f'p_value_{type_test}'] < row[f'(i/m)a_{type_test}']) and row[f'rank_{type_test}'] == (rank_begin+1):
            rank_begin = row[f'rank_{type_test}']
            list_true_false.append(True)
        else:
            list_true_false.append(False)
    df[f'check_bh_{type_test}'] = list_true_false
    return df
    


def get_significant(df, type_test):
    df_select = df[['info', f'p_value_{type_test}', f'bonferroni_{type_test}', f'fdr_bh_{type_test}']]
    print(f'---------{type_test}--------')
    normal_p = df_select[df_select[f'p_value_{type_test}'] == True]
    genes_sig = normal_p['info']
    print(f'nomal: {len(genes_sig)}')

    bon_p = df_select[df_select[f'bonferroni_{type_test}'] == True]
    bon_sig = bon_p['info']
    print(f'bonferroni: {len(bon_sig)}')

    bh_p = df_select[df_select[f'fdr_bh_{type_test}'] == True]
    bh_sig = bh_p['info']
    print(f'fdr_bh: {len(bh_sig)}')
    return genes_sig, bon_sig, bh_sig


def get_overlap(self_X2_df, X2_df, F_df):
    print('---ALL')
    all_values = set(self_X2_df) | set(X2_df) |  set(F_df)
    # print(f'All genes ({len(all_values)}): {all_values}')
    print(f'All genes ({len(all_values)})')

    elements_in_all = list(set.intersection(*map(set, [self_X2_df, X2_df, F_df])))
    # print(f'Gene in all list ({len(elements_in_all)}): {elements_in_all}\n')
    print(f'Gene in all list ({len(elements_in_all)})')

    # selfX2 and X2
    print('---selfX2 and X2')
    in_selfX2_X2 = set(self_X2_df) & set(X2_df)
    # print(f'in selfX2 and X2 ({len(in_selfX2_X2)}): {in_selfX2_X2}')
    print(f'in selfX2 and X2 ({len(in_selfX2_X2)})')

    dif_selfX2_X2 = set(self_X2_df) - set(X2_df)
    # print(f'in selfX2 NOT in X2 ({len(dif_selfX2_X2)}): {dif_selfX2_X2}')
    print(f'in selfX2 NOT in X2 ({len(dif_selfX2_X2)})')

    dif_X2_selfX2 = set(X2_df) - set(self_X2_df)
    # print(f'in X2 NOT in selfx2 ({len(dif_X2_selfX2)}): {dif_X2_selfX2}\n')
    print(f'in X2 NOT in selfx2 ({len(dif_X2_selfX2)})')

    # selfX2 and F
    print('---selfX2 and F')
    in_selfX2_F = set(self_X2_df) & set(F_df)
    # print(f'in selfX2 and F ({len(in_selfX2_F)}): {in_selfX2_F}')
    print(f'in selfX2 and F ({len(in_selfX2_F)})')

    dif_selfX2_F = set(self_X2_df) - set(F_df)
    # print(f'in selfX2 NOT in F ({len(dif_selfX2_F)}): {dif_selfX2_F}')
    print(f'in selfX2 NOT in F ({len(dif_selfX2_F)})')

    dif_F_selfX2 = set(F_df) - set(self_X2_df)
    # print(f'in F NOT in selfX2 ({len(dif_F_selfX2)}): {dif_F_selfX2}\n')
    print(f'in F NOT in selfX2 ({len(dif_F_selfX2)})')

    # F and X2
    print('---F and X2')
    in_F_X2 = set(F_df) & set(X2_df)
    # print(f'in F and X2 ({len(in_F_X2)}): {in_F_X2}')
    print(f'in F and X2 ({len(in_F_X2)})')

    dif_F_X2 = set(F_df) - set(X2_df)
    # print(f'in F NOT in X2 ({len(dif_F_X2)}): {dif_F_X2}')
    print(f'in F NOT in X2 ({len(dif_F_X2)})')

    dif_X2_F = set(X2_df) - set(F_df)
    # print(f'in X2 NOT in F ({len(dif_X2_F)}): {dif_X2_F}')
    print(f'in X2 NOT in F ({len(dif_X2_F)})\n')
    return elements_in_all

def search(df_select, type_file, non_coding, path_analyse, GT, fc):
    alpha=0.05
    if GT:
        method = 'bonferroni'
        df_select = bonferroni_correction(df_select, alpha, method, 'cochran_armitage')
        method = 'fdr_bh'
        df_select = benjamini_hochberg_correction(df_select, alpha, 'cochran_armitage')
        
        # df_select.to_csv(f"{path_analyse}{type_file}_{non_coding}_{fc}_MTC.tsv", sep='\t', encoding='utf-8', index=False)

        cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig = get_significant(df_select, 'cochran_armitage', alpha)
        elements_snps_all_MTC = list(set.intersection(*map(set, [cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig])))
        return cochran_armitage_sig_normal, cochran_armitage_bon_sig, cochran_armitage_bh_sig, elements_snps_all_MTC
    else:
        method = 'bonferroni'
        df_select = bonferroni_correction(df_select, alpha, method, 'X2_self')
        df_select = bonferroni_correction(df_select, alpha, method, 'X2')
        df_select = bonferroni_correction(df_select, alpha, method, 'F')
        
        method = 'fdr_bh'
        df_select = benjamini_hochberg_correction(df_select, alpha, 'X2_self')
        df_select = benjamini_hochberg_correction(df_select, alpha, 'X2')
        df_select = benjamini_hochberg_correction(df_select, alpha, 'F')
        print(df_select.head())
        print(df_select.shape)
        print(len(df_select))

        # df_select.to_csv(f"{path_analyse}{type_file}_{non_coding}_{fc}_MTC.tsv", sep='\t', encoding='utf-8', index=False)

        X2_self_sig_normal, X2_self_bon_sig, X2_self_bh_sig = get_significant(df_select, 'X2_self', alpha)
        X2_sig_normal, X2_bon_sig, X2_bh_sig = get_significant(df_select, 'X2', alpha)
        F_sig_normal, F_bon_sig, F_bh_sig = get_significant(df_select, 'F', alpha)

        print('\nNORMAL')
        elements_in_all_normal = get_overlap(X2_self_sig_normal, X2_sig_normal, F_sig_normal) #[:top_num]
        print('\nBON')
        elements_in_all_bon = get_overlap(X2_self_bon_sig, X2_bon_sig, F_bon_sig)
        print('\nBH')
        elements_in_all_bh = get_overlap(X2_self_bh_sig, X2_bh_sig, F_bh_sig)
        elements_snps_all_MTC = list(set.intersection(*map(set, [elements_in_all_normal, elements_in_all_bon, elements_in_all_bh])))
        print(elements_snps_all_MTC)
        return elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC