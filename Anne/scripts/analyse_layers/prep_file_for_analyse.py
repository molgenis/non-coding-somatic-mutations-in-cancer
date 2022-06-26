#!/usr/bin/env python3

# Imports
import sys
import pandas as pd

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
import ast


def prep_file(path_file, path_save, type_data):
    """
    Filter the file and make it so that it can be used in the analyses.
    :param path_file:  Path to file
    :param path_save:  Path where the new file should be saved
    :param type_data:  Type of data/name of layer
    :return: df_b_nb: Dataframe with breast cancer counts and non breast cancer counts
    """
    df = pd.read_csv(path_file, sep='\t')
    # Make list
    breast_count_list = list()
    nonbreast_count_list = list()
    # Loop over rows
    for index, row in df.iterrows():
        # Keep track of how many counts there are for breast cancer and not breast cancer
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
    df_b_nb.to_csv(f"{path_save}b_nb_gene_{type_data}_2000_250.tsv", sep='\t', encoding='utf-8', index=False)
    return df_b_nb


def prep_file_multi(path_file, path_save, type_data, select_chrom):
    """
    Filter the file and make it so that it can be used in the analyses (multiprocess).
    :param path_file:  Path to file
    :param path_save:  Path where the new file should be saved
    :param type_data:  Type of data/name of layer
    :param select_chrom: The selected chromosome 
    :return: df_b_nb: Dataframe with breast cancer counts and non breast cancer counts
    """
    df = pd.read_csv(path_file, sep='\t')
    # Filter df on chromosome
    df = df[df['chr'] == select_chrom]
    # Make list
    breast_count_list = list()
    nonbreast_count_list = list()
    # Loop over rows
    for index, row in df.iterrows():
        # Keep track of how many counts there are for breast cancer and not breast cancer
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
    df_b_nb.to_csv(f"{path_save}b_nb_{type_data}_{select_chrom}.tsv", sep='\t', encoding='utf-8', index=False)
    return df_b_nb
