import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import search_close_gene as search_close_gene
import bon_and_bh_calculate as bon_and_bh_calculate
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



def get_significant_snps(list_significant_elements, f, type_test):
    print(len(list_significant_elements))
    significant_snps = list()
    for i in list_significant_elements:
        # Append snp_ID
        significant_snps.append(i.split('_')[0])
    f.write(f"{type_test}\t{','.join(map(str, list(significant_snps)))}\n")
    return significant_snps



def run_different_fc(df_select, type_file, non_coding, path_analyse, fc):
    elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC = bon_and_bh_calculate.search(df_select, type_file, non_coding, path_analyse, True, fc)
    
    search_close_gene.search_gene(elements_in_all_normal, path_analyse, type_file, non_coding, 'normal', fc)
    search_close_gene.search_gene(elements_in_all_bon, path_analyse, type_file, non_coding, 'bon', fc)
    search_close_gene.search_gene(elements_in_all_bh, path_analyse, type_file, non_coding, 'bh', fc)
    search_close_gene.search_gene(elements_snps_all_MTC, path_analyse, type_file, non_coding, 'all_MTC', fc)

    
    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{fc}_sig_snps.tsv", "a")
    significant_snps_normal = get_significant_snps(elements_in_all_normal, f, 'snps_normal')
    significant_snps_bon = get_significant_snps(elements_in_all_bon, f, 'snps_bon')
    significant_snps_bh = get_significant_snps(elements_in_all_bh, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps(elements_snps_all_MTC, f, 'snps_all_MTC')
    f.close()



def run_all_corrections(path_analyse, type_file, non_coding):
    print(type_file)
    print(non_coding)
    path_file = f"{path_analyse}{type_file}_{non_coding}_cochran_armitage.tsv"
    df = pd.read_csv(path_file, sep='\t')
    print(df.columns)
    df_select = df[['snp_ID', 'chr', 'pos_start', 'pos_end', 'GT_1_b', 'GT_2_b', 'GT_0_b', 'GT_1_nb', 'GT_2_nb', 'GT_0_nb', 'p_value_cochran_armitage']]
    df_select['foldchange'] = (df_select['GT_1_b']+(df_select['GT_2_b']*2))/(df_select['GT_1_nb']+(df_select['GT_2_nb']*2))
    df_select['info'] = df_select['snp_ID'].map(str) + '_' + df_select['chr'].map(str) + '_' + df_select['pos_start'].map(str) + '_' + df_select['pos_end'].map(str) + '_' + df_select['foldchange'].map(str)

    run_different_fc(df_select, type_file, non_coding, path_analyse, 'ALL')

    df_breast = df_select[df_select['foldchange'] > 1]
    run_different_fc(df_breast, type_file, non_coding, path_analyse, 'breast')

    df_nonbreast = df_select[df_select['foldchange'] <= 1]
    run_different_fc(df_nonbreast, type_file, non_coding, path_analyse, 'nonbreast')


   
    

def main():
    config = get_config('gearshift')
    path_analyse = config['analyse_new']

    # per_snp
    type_file = 'GT'
    non_coding = 'ALL'
    run_all_corrections(path_analyse, type_file, non_coding)
    non_coding = 'Coding'
    run_all_corrections(path_analyse, type_file, non_coding)
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding)



if __name__ == '__main__':
    main()