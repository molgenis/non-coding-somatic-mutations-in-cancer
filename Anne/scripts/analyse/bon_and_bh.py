import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests


def run_multiple_testing_correction(df, alpha, method, type_test):
    list_p_values = list(df[f'p_value_{type_test}'])
    rej, p_adjusted = multipletests(list_p_values, alpha=alpha, method=method)[:2]
    df[f'{method}_{type_test}'] = p_adjusted
    return df

def get_significant(df, type_test, alpha):
    df_select = df[['info', f'p_value_{type_test}', f'bonferroni_{type_test}', f'fdr_bh_{type_test}']]
    print(f'---------{type_test}--------')
    normal_p = df_select[df_select[f'p_value_{type_test}'] <= alpha]
    sort_normal_p = normal_p.sort_values(f'p_value_{type_test}')
    genes_sig = sort_normal_p['info']
    print(f'nomal: {len(genes_sig)}')

    bon_p = df_select[df_select[f'bonferroni_{type_test}'] <= alpha]
    sort_bon = bon_p.sort_values(f'bonferroni_{type_test}')
    bon_sig = sort_bon['info']
    print(f'bonferroni: {len(bon_sig)}')

    bh_p = df_select[df_select[f'fdr_bh_{type_test}'] <= alpha]
    sort_bh = bh_p.sort_values(f'fdr_bh_{type_test}')
    bh_sig = sort_bh['info']
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
        df_select['info'] = df_select['chr'].map(str) + '__' + df_select['start_position_regio'].map(str) + '__' + df_select['end_position_regio'].map(str)
    print(df_select)

    alpha=0.05
    method = 'bonferroni'
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'X2_self')
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'X2')
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'F')
    
    method = 'fdr_bh'
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'X2_self')
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'X2')
    df_select = run_multiple_testing_correction(df_select, alpha, method, 'F')
    print(df_select)
    print(len(df_select))

    # df_select.to_csv(f"{path_analyse}{type_file}_{non_coding}_MTC.tsv", sep='\t', encoding='utf-8', index=False)

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

    # path_search_snp = "D:/Hanze_Groningen/STAGE/analyse/stat/ALL_gene_before_2000_250.tsv"
    snps_search = pd.read_csv(path_search_snp, sep='\t')
    if with_gene:
        snps_search['info'] = snps_search['gene'].map(str) + '__' + snps_search['chr'].map(str) + '__' + snps_search['start_position_regio'].map(str) + '__' + snps_search['end_position_regio'].map(str)
    else:
        snps_search['info'] = snps_search['chr'].map(str) + '__' + snps_search['start_position_regio'].map(str) + '__' + snps_search['end_position_regio'].map(str)
    
    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_sig_snps.tsv", "a")
    significant_snps_normal = get_significant_snps(elements_in_all_normal, snps_search, f, 'snps_normal')
    significant_snps_bon = get_significant_snps(elements_in_all_bon, snps_search, f, 'snps_bon')
    significant_snps_bh = get_significant_snps(elements_in_all_bh, snps_search, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps(elements_snps_all_MTC, snps_search, f, 'snps_all_MTC')
    f.close()
    

def main():
    path_analyse = 'D:/Hanze_Groningen/STAGE/analyse/new/'
    path_snp_ids = 'D:/Hanze_Groningen/STAGE/lagen/'
    with_gene = False
    # TFBS, UCNE, DNase
    type_file = 'TFBS'
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    # # type_file = 'UCNE'
    # # run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    # # type_file = 'DNase'
    # # run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}{type_file}_chrALL_num_snps.tsv')
    
    # AfterGene/BeforeGene
    with_gene = True
    type_file = 'afterGene'
    non_coding = 'NonCoding_Coding'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}ALL_gene_after_2000_250.tsv')
    type_file = 'beforeGene'
    run_all_corrections(path_analyse, type_file, non_coding, with_gene, f'{path_snp_ids}ALL_gene_before_2000_250.tsv')






if __name__ == '__main__':
    main()