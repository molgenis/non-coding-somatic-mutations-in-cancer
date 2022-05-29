import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
from Database import Database


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

def get_significant_snps(list_significant_elements, db, f, type_test):
    print(len(list_significant_elements))
    snp_id_list = list()
    dict_snp = dict()
    for sig_ele in list_significant_elements:
        print(sig_ele)
        db.cursor.execute("""
                        SELECT ID
                        FROM 'snp'
                        WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                        """ %
                        (str(sig_ele.split('_')[0].replace('chr', '')), int(sig_ele.split('_')[1]), int(sig_ele.split('_')[2])))
        results = db.cursor.fetchall()
        # Make snp_id_list
        list_element = list()
        for res in results:
            # Add ID to snp_id_list
            snp_id_list.append(res['ID'])
            list_element.append(res['ID'])
        dict_snp[sig_ele] = list_element

    f.write(f"{type_test}\t{','.join(map(str, list(set(snp_id_list))))}\t{dict_snp}\n")



def run_all_corrections(path_analyse, type_file, non_coding, db, col1, col2):
    print(type_file)
    print(non_coding)
    path_file = f"{path_analyse}{type_file}_{non_coding}_both_0_TESTS.tsv"
    df = pd.read_csv(path_file, sep='\t')
    print(df.columns)
    df_select = df[['chr', col1, col2, 'counts_breast', 'counts_nonbreast',
                'p_value_X2_self', 'p_value_X2', 'p_value_F']]
    df_select['info'] = df_select['chr'].map(str) + '_' + df_select[col1].map(str) + '_' + df_select[col2].map(str)

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

    # df_select.to_csv(f"{path_analyse}correction/{type_file}_{non_coding}_MTC.tsv", sep='\t', encoding='utf-8', index=False)

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

   
    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_sig_snps.tsv", "a")
    significant_snps_normal = get_significant_snps(elements_in_all_normal, db, f, 'snps_normal')
    significant_snps_bon = get_significant_snps(elements_in_all_bon, db, f, 'snps_bon')
    significant_snps_bh = get_significant_snps(elements_in_all_bh, db, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps(elements_snps_all_MTC, db, f, 'snps_all_MTC')
    f.close()
    

def main():
    path_db = 'D:/Hanze_Groningen/STAGE/lastdb/db_laatste_copy.db' #config['database']
    db = Database(path_db)
    path_analyse = 'D:/Hanze_Groningen/STAGE/analyse/new/'

    # per_snp
    type_file = 'per_snp'
    non_coding = 'ALL'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'pos_start', 'pos_end')
    non_coding = 'Coding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'pos_start', 'pos_end')
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'pos_start', 'pos_end')

    # Region 1000
    type_file = 'Region_1000'
    non_coding = 'ALL'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')
    non_coding = 'Coding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')

    # Region 2000
    type_file = 'Region_2000'
    non_coding = 'ALL'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')
    non_coding = 'Coding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')
    non_coding = 'NonCoding'
    run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')





if __name__ == '__main__':
    main()