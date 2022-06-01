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



def get_significant_snps_db(list_significant_elements, db, f, type_test):
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
                        (str(sig_ele.split('_')[0].replace('chr', '')), int(round(float(sig_ele.split('_')[1]))), int(round(float(sig_ele.split('_')[2])))))
                        
        results = db.cursor.fetchall()
        # Make snp_id_list
        list_element = list()
        for res in results:
            # Add ID to snp_id_list
            snp_id_list.append(res['ID'])
            list_element.append(res['ID'])
        dict_snp[sig_ele] = list_element

    f.write(f"{type_test}\t{','.join(map(str, list(set(snp_id_list))))}\t{dict_snp}\n")


def run_different_fc(df_select, type_file, non_coding, path_analyse, db, fc):
    elements_in_all_normal, elements_in_all_bon, elements_in_all_bh, elements_snps_all_MTC = bon_and_bh_calculate.search(df_select, type_file, non_coding, path_analyse, False, fc)
    

    search_close_gene.search_gene(elements_in_all_normal, path_analyse, type_file, non_coding, 'normal', fc)
    search_close_gene.search_gene(elements_in_all_bon, path_analyse, type_file, non_coding, 'bon', fc)
    search_close_gene.search_gene(elements_in_all_bh, path_analyse, type_file, non_coding, 'bh', fc)
    search_close_gene.search_gene(elements_snps_all_MTC, path_analyse, type_file, non_coding, 'all_MTC', fc)

   
    # Open and create file
    f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{fc}_sig_snps.tsv", "a")
    significant_snps_normal = get_significant_snps_db(elements_in_all_normal, db, f, 'snps_normal')
    significant_snps_bon = get_significant_snps_db(elements_in_all_bon, db, f, 'snps_bon')
    significant_snps_bh = get_significant_snps_db(elements_in_all_bh, db, f, 'snps_bh')
    significant_snps_all_MTC = get_significant_snps_db(elements_snps_all_MTC, db, f, 'snps_all_MTC')
    f.close()


def run_all_corrections(path_analyse, type_file, non_coding, db, col1, col2):
    print(type_file)
    print(non_coding)
    path_file = f"{path_analyse}{type_file}_{non_coding}_both_0_TESTS.tsv"
    df = pd.read_csv(path_file, sep='\t')
    print(df.columns)
    df_select = df[['chr', col1, col2, 'counts_breast', 'counts_nonbreast',
                'p_value_X2_self', 'p_value_X2', 'p_value_F']]
    df_select['foldchange'] = df_select['counts_breast']/df_select['counts_nonbreast']
    df_select['info'] = df_select['chr'].map(str) + '_' + df_select[col1].map(str) + '_' + df_select[col2].map(str) + '_' + df_select['foldchange'].map(str)

    run_different_fc(df_select, type_file, non_coding, path_analyse, db, 'ALL')

    df_breast = df_select[df_select['foldchange'] > 1]
    run_different_fc(df_breast, type_file, non_coding, path_analyse, db, 'breast')

    df_nonbreast = df_select[df_select['foldchange'] <= 1]
    run_different_fc(df_nonbreast, type_file, non_coding, path_analyse, db, 'nonbreast')


def main():
    config = get_config('gearshift')
    path_db = config['database']
    db = Database(path_db)
    path_analyse = config['analyse_new']

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
    db.close()





if __name__ == '__main__':
    main()