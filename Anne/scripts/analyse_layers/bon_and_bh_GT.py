#!/usr/bin/env python3

#Imports
import pandas as pd
import search_close_gene as search_close_gene
import bon_and_bh_calculate as bon_and_bh_calculate
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def get_significant_snps(list_significant_elements, f, type_test):
    """
    Retrieves the SNP IDs from the significant regions
    :param list_significant_elements: List of significant regions
    :param f:   File in which line is written
    :param type_test:    Which test manual chi-square, chi-square or fisher exactly test    
    :return:    significant_snps: List of SNP ids
    """
    significant_snps = list()
    for i in list_significant_elements:
        # Append snp_ID
        significant_snps.append(i.split('_')[0])
    f.write(f"{type_test}\t{','.join(map(str, list(significant_snps)))}\n")
    return significant_snps


def run_different_fc(df_select, type_file, non_coding, path_analyse, fc):
    """
    Calls functions. These functions search for the significant regions, the closest genes to these regions
    and the SNP IDs within these regions
    :param df_select:  The dataframe with the p-values
    :param type_file:  What type it is, TFBS, GT or any of the other options
    :param non_coding:  Whether they are non-coding areas or coding or a mix of both
    :param path_analyse:  The path where the new files will be saved
    :param fc:      Which foldchange is used        
    :return:    
    """
    # Makes the multiple testing corrections and extracts the significant regions here
    elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC = bon_and_bh_calculate.search(
        df_select, type_file, non_coding, path_analyse, True, fc)
    # Finds the closest gene in the significant regions
    search_close_gene.search_gene(elements_in_all_normal, path_analyse, type_file, non_coding, 'normal', fc)
    search_close_gene.search_gene(elements_in_all_bon, path_analyse, type_file, non_coding, 'bon', fc)
    search_close_gene.search_gene(elements_in_all_bh, path_analyse, type_file, non_coding, 'bh', fc)
    search_close_gene.search_gene(elements_snps_all_MTC, path_analyse, type_file, non_coding, 'all_MTC', fc)

    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{fc}_sig_snps.tsv", "a")
    # Grabs all SNP IDs that fall within the significant regions
    significant_snps_normal = get_significant_snps(elements_in_all_normal, f, 'snps_normal')
    significant_snps_bon = get_significant_snps(elements_in_all_bon, f, 'snps_bon')
    significant_snps_bh = get_significant_snps(elements_in_all_bh, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps(elements_snps_all_MTC, f, 'snps_all_MTC')
    f.close()


def run_all_corrections(path_analyse, type_file, non_coding):
    """
    Calculates the foldchange and calls run_different_fc
    :param path_analyse:  The path where the new files will be saved    
    :param type_file:     What type it is, TFBS, GT or any of the other options
    :param non_coding:    Whether they are non-coding areas or coding or a mix of both
    :return:    
    """
    # Read file
    path_file = f"{path_analyse}{type_file}_{non_coding}_cochran_armitage.tsv.gz"
    df = pd.read_csv(path_file, sep='\t', compression='gzip')
    # Select df
    df_select = df[
        ['snp_ID', 'chr', 'pos_start', 'pos_end', 'GT_1_b', 'GT_2_b', 'GT_0_b', 'GT_1_nb', 'GT_2_nb', 'GT_0_nb',
         'p_value_cochran_armitage']]
    # Calculate foldchange
    df_select['foldchange'] = (df_select['GT_1_b'] + (df_select['GT_2_b'] * 2)) / (
                df_select['GT_1_nb'] + (df_select['GT_2_nb'] * 2))
    df_select['info'] = df_select['snp_ID'].map(str) + '_' + df_select['chr'].map(str) + '_' + df_select[
        'pos_start'].map(str) + '_' + df_select['pos_end'].map(str) + '_' + df_select['foldchange'].map(str) + '_None'

    run_different_fc(df_select, type_file, non_coding, path_analyse, 'ALL')


def main():
    # Call get_config
    config = get_config('gearshift')
    path_analyse = config['analyse_new']

    # per_snp
    type_file = 'GT'
    for non_coding in ['ALL', 'Coding', 'NonCoding']:
        run_all_corrections(path_analyse, type_file, non_coding)


if __name__ == '__main__':
    main()
