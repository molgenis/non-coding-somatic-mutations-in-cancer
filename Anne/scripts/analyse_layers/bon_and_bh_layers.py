#!/usr/bin/env python3

#Imports
import pandas as pd
import numpy as np
import search_close_gene as search_close_gene
import bon_and_bh_calculate as bon_and_bh_calculate
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def get_significant_snps(df, list_significant_elements, snps_search, path_analyse, type_test, type_file, non_coding,
                         fc):
    """
    Retrieves the SNP IDs from the significant regions
    :param df: The dataframe
    :param list_significant_elements:  List of significant regions
    :param snps_search: File from which the SNP ids can be obtained
    :param path_analyse: The path where the new files will be saved
    :param type_test:  Which test manual chi-square, chi-square or fisher exactly test  
    :param type_file: What type it is, TFBS, GT or any of the other options
    :param non_coding: Whether they are non-coding areas or coding or a mix of both
    :param fc:  Which foldchange is used             
    :return:    significant_snps: List of SNP ids
    """
    path_new_file = f"{path_analyse}correction/{type_file}_{non_coding}_{fc}_{type_test}_sig_snps.tsv"
    snp_DF = pd.DataFrame(columns=[f'snp_IDs_{type_test}'])
    # Check if significant regions were found  
    if len(list_significant_elements) > 0:
        info_chr_list = list(df[df['info'].isin(list_significant_elements)]['info_chr'])
        significant_elements_df = snps_search[snps_search['info_chr'].isin(info_chr_list)]
        significant_snps = list(set(','.join(list(significant_elements_df['snp_list'])).split(',')))
        snp_DF[f'snp_IDs_{type_test}'] = significant_snps
        snp_DF.to_csv(path_new_file, sep='\t', encoding='utf-8', index=False)
        return significant_snps
    else:
        snp_DF[f'snp_IDs_{type_test}'] = ['NONE']
        snp_DF.to_csv(path_new_file, sep='\t', encoding='utf-8', index=False)
        return []


def run_different_fc(df_select, type_file, non_coding, path_analyse, fc, path_search_snp, with_gene):
    """
    Calls functions. These functions search for the significant regions, the closest genes to these regions
    and the SNP IDs within these regions
    :param df_select:  The dataframe with the p-values
    :param type_file:  What type it is, TFBS, GT or any of the other options
    :param non_coding:  Whether they are non-coding areas or coding or a mix of both
    :param path_analyse:  The path where the new files will be saved
    :param fc:      Which foldchange is used  
    :param path_search_snp: Path to file from which the SNP ids can be obtained
    :param with_gene:        Whether it's before/after a gene
    :return:    
    """
    # Makes the multiple testing corrections and extracts the significant regions here
    elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC = bon_and_bh_calculate.search(
        df_select, type_file, non_coding, path_analyse, False, fc)

    snps_search = pd.read_csv(path_search_snp, sep='\t', compression='gzip')
    if with_gene:
        snps_search['info_chr'] = snps_search['gene'].map(str) + '__' + snps_search['chr'].map(str) + '__' + \
                                  snps_search['start_position_regio'].map(str) + '__' + snps_search[
                                      'end_position_regio'].map(str)
    else:
        snps_search['info_chr'] = snps_search['chr'].map(str) + '_' + snps_search['start_position_regio'].map(
            str) + '_' + snps_search['end_position_regio'].map(str)
    # Finds the closest gene in the significant regions
    search_close_gene.search_gene(elements_in_all_normal, path_analyse, type_file, non_coding, 'normal', fc)
    search_close_gene.search_gene(elements_in_all_bon, path_analyse, type_file, non_coding, 'bon', fc)
    search_close_gene.search_gene(elements_in_all_bh, path_analyse, type_file, non_coding, 'bh', fc)
    search_close_gene.search_gene(elements_snps_all_MTC, path_analyse, type_file, non_coding, 'all_MTC', fc)

    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{fc}_sig_snps.tsv", "a")
    # Grabs all SNP IDs that fall within the significant regions
    significant_snps_normal = get_significant_snps(df_select, elements_in_all_normal, snps_search, path_analyse,
                                                   'snps_normal', type_file, non_coding, fc)
    significant_snps_bon = get_significant_snps(df_select, elements_in_all_bon, snps_search, path_analyse, 'snps_bon',
                                                type_file, non_coding, fc)
    significant_snps_bh = get_significant_snps(df_select, elements_in_all_bh, snps_search, path_analyse, 'snps_bh',
                                               type_file, non_coding, fc)
    significant_snps_all_MTC = get_significant_snps(df_select, elements_snps_all_MTC, snps_search, path_analyse,
                                                    'snps_all_MTC', type_file, non_coding, fc)
    f.close()


def calculate_fc(df):
    """
    Calculate the fold change
    :param df: The dataframe       
    :return: df: with foldchange column   
    """
    df['foldchange'] = df['counts_breast'] / df['counts_nonbreast']
    # When you divide 0 by 0 you get a nan value. We therefore replace these with 1, because they
    # are equal to each other.
    df['bigger'] = df['foldchange'].fillna(1)
    # When you divide a number by 0 you get an 'inf' value. This is replaced with 2 because the breast
    # cancer counts are higher.
    df['bigger'].replace([np.inf, -np.inf], 2, inplace=True)
    change_bigger = np.array(df['bigger'].values.tolist())
    df['bigger'] = np.where(change_bigger > 1.0, 2, change_bigger).tolist()
    change_bigger = np.array(df['bigger'].values.tolist())
    # When you divide 0 by a number, you get 0. In this case, the non-breast cancer is larger.
    df['bigger'] = np.where(change_bigger < 1.0, 0, change_bigger).tolist()
    df['bigger'].replace(0, 'nb', inplace=True)
    df['bigger'].replace(1, '=', inplace=True)
    df['bigger'].replace(2, 'b', inplace=True)
    return df


def run_all_corrections(path_analyse, type_file, non_coding, with_gene, path_search_snp):
    """
    Calculates the foldchange and calls run_different_fc
    :param path_analyse:  The path where the new files will be saved    
    :param type_file:  What type it is, TFBS, GT or any of the other options
    :param non_coding:  Whether they are non-coding areas or coding or a mix of both 
    :param with_gene:  Whether it's before/after a gene
    :param path_search_snp: Path to file from which the SNP ids can be obtained 
    :return:    
    """
    path_file = f"{path_analyse}{type_file}_{non_coding}_both_0_TESTS_chr0_0.tsv.gz"
    df = pd.read_csv(path_file, sep='\t', compression='gzip')

    if with_gene:
        df_select = df[
            ['gene', 'chr', 'start_position_regio', 'end_position_regio', 'counts_breast', 'counts_nonbreast',
             'p_value_X2_self', 'p_value_X2', 'p_value_F', 'one_in_interval_bon']]
        df_select = calculate_fc(df_select)
        df_select['info'] = df_select['gene'].map(str) + '__' + df_select['chr'].map(str) + '__' + df_select[
            'start_position_regio'].map(str) + '__' + df_select['end_position_regio'].map(str) + '__' + df_select[
                                'foldchange'].map(str) + '__' + df_select['bigger'].map(str)
        df_select['info_chr'] = df_select['gene'].map(str) + '__' + df_select['chr'].map(str) + '__' + df_select[
            'start_position_regio'].map(str) + '__' + df_select['end_position_regio'].map(str)

    else:
        df_select = df[['chr', 'start_position_regio', 'end_position_regio', 'counts_breast', 'counts_nonbreast',
                        'p_value_X2_self', 'p_value_X2', 'p_value_F', 'one_in_interval_bon']]
        df_select = calculate_fc(df_select)
        df_select['info'] = df_select['chr'].map(str) + '_' + df_select['start_position_regio'].map(str) + '_' + \
                            df_select['end_position_regio'].map(str) + '_' + df_select['foldchange'].map(str) + '_' + \
                            df_select['bigger'].map(str)
        df_select['info_chr'] = df_select['chr'].map(str) + '_' + df_select['start_position_regio'].map(str) + '_' + \
                                df_select['end_position_regio'].map(str)

    # For all significant regions
    run_different_fc(df_select, type_file, non_coding, path_analyse, 'ALL', path_search_snp, with_gene)
    # For all significant regions with a fc > 1 
    df_breast = df_select[df_select['bigger'] == 'b']
    run_different_fc(df_breast, type_file, non_coding, path_analyse, 'breast', path_search_snp, with_gene)
    # For all significant regions with a fc < 1
    df_nonbreast = df_select[df_select['bigger'] == 'nb']
    run_different_fc(df_nonbreast, type_file, non_coding, path_analyse, 'nonbreast', path_search_snp, with_gene)


def main():
    # Call get_config
    config = get_config('gearshift')
    path_analyse = config['analyse_new']
    path_snp_ids = config['layers']
    # Loop over layers
    for type_file in ['TFBS', 'UCNE', 'DNase', 'afterGene', 'beforeGene']:
        if 'Gene' in type_file:
            # AfterGene/BeforeGene
            with_gene = True
            non_coding = 'NonCoding_Coding'
            run_all_corrections(path_analyse, type_file, non_coding, with_gene,
                                f"{path_snp_ids}ALL_gene_{type_file.replace('Gene', '')}_2000_250.tsv.gz")
        else:
            # TFBS, UCNE, DNase
            with_gene = False
            non_coding = 'NonCoding2'
            run_all_corrections(path_analyse, type_file, non_coding, with_gene,
                                f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv.gz')
            non_coding = 'NonCoding_NC2'
            run_all_corrections(path_analyse, type_file, non_coding, with_gene,
                                f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv.gz')


if __name__ == '__main__':
    main()
