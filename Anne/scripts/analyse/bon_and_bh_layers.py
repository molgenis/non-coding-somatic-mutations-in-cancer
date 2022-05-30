import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import search_close_gene as search_close_gene
import bon_and_bh_calculate as bon_and_bh_calculate
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



def get_significant_snps(list_significant_elements, df, f, type_test):
    print(len(list_significant_elements))
    significant_elements_df = df[df['info'].isin(list_significant_elements)]
    significant_snps = list(set(','.join(list(significant_elements_df['snp_list'])).split(',')))
    print(significant_snps)
    print(len(significant_snps))
    print()
    f.write(f"{type_test}\t{','.join(map(str, list(significant_snps)))}\n")
    return significant_snps



def run_all_corrections(path_analyse, type_file, non_coding, with_gene, path_search_snp):
    path_file = f"{path_analyse}{type_file}_{non_coding}_both_0_TESTS.tsv"
    df = pd.read_csv(path_file, sep='\t')
    print(df.head())
    if with_gene:
        df_select = df[['gene', 'chr', 'start_position_regio', 'end_position_regio', 'counts_breast', 'counts_nonbreast',
                    'p_value_X2_self', 'p_value_X2', 'p_value_F']]
        df_select['info'] = df_select['gene'].map(str) + '__' + df_select['chr'].map(str) + '__' + df_select['start_position_regio'].map(str) + '__' + df_select['end_position_regio'].map(str)
    else:
        df_select = df[['chr', 'start_position_regio', 'end_position_regio', 'counts_breast', 'counts_nonbreast',
                    'p_value_X2_self', 'p_value_X2', 'p_value_F']]
        df_select['info'] = df_select['chr'].map(str) + '_' + df_select['start_position_regio'].map(str) + '_' + df_select['end_position_regio'].map(str)
    print(df_select)

    elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC = bon_and_bh_calculate.search(df_select, type_file, non_coding, path_analyse, False)

    # path_search_snp = "D:/Hanze_Groningen/STAGE/analyse/stat/ALL_gene_before_2000_250.tsv"
    snps_search = pd.read_csv(path_search_snp, sep='\t')
    if with_gene:
        snps_search['info'] = snps_search['gene'].map(str) + '__' + snps_search['chr'].map(str) + '__' + snps_search['start_position_regio'].map(str) + '__' + snps_search['end_position_regio'].map(str)
    else:
        snps_search['info'] = snps_search['chr'].map(str) + '_' + snps_search['start_position_regio'].map(str) + '_' + snps_search['end_position_regio'].map(str)
        search_close_gene.search_gene(elements_in_all_normal, path_analyse, type_file, non_coding, 'normal')
        search_close_gene.search_gene(elements_in_all_bon, path_analyse, type_file, non_coding, 'bon')
        search_close_gene.search_gene(elements_in_all_bh, path_analyse, type_file, non_coding, 'bh')
        search_close_gene.search_gene(elements_snps_all_MTC, path_analyse, type_file, non_coding, 'all_MTC')

    
    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_sig_snps.tsv", "a")
    significant_snps_normal = get_significant_snps(elements_in_all_normal, snps_search, f, 'snps_normal')
    significant_snps_bon = get_significant_snps(elements_in_all_bon, snps_search, f, 'snps_bon')
    significant_snps_bh = get_significant_snps(elements_in_all_bh, snps_search, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps(elements_snps_all_MTC, snps_search, f, 'snps_all_MTC')
    f.close()
    

def main():
    config = get_config()
    path_analyse = config['analyse_new'] #'D:/Hanze_Groningen/STAGE/analyse/new/' #config['analyse_new']
    path_snp_ids = config['layers'] #config['layers'] #'D:/Hanze_Groningen/STAGE/lagen/'
    with_gene = False
    # TFBS, UCNE, DNase
    type_file = 'TFBS'
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    type_file = 'UCNE'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    type_file = 'DNase'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    
    # AfterGene/BeforeGene
    with_gene = True
    type_file = 'afterGene'
    non_coding = 'NonCoding_Coding'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}ALL_gene_after_2000_250.tsv')
    type_file = 'beforeGene'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}ALL_gene_before_2000_250.tsv')






if __name__ == '__main__':
    main()