import pandas as pd
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import search_close_gene as search_close_gene
import bon_and_bh_calculate as bon_and_bh_calculate
import sys
# sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
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


def calculate_fc(df):
    """

    example:
    np.array([10, 20, 0, 0, 3, 4]) / np.array([0, 20, 0, 30, 2, 33])   = array([inf, 1. , nan, 0. , 1.5, 0.12121212])

    """
    df['foldchange'] = df['counts_breast']/df['counts_nonbreast']
    # When you divide 0 by 0 you get a nan value. We therefore replace these with 1, because they are equal to each other.
    df['bigger'] = df['foldchange'].fillna(1)
    # When you divide a number by 0 you get an 'inf' value. This is replaced with 2 because the breast cancer counts are higher.
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


def run_all_corrections(path_analyse, type_file, non_coding, db, col1, col2):
    print(type_file)
    print(non_coding)
    path_file = f"{path_analyse}{type_file}_{non_coding}_both_0_TESTS_chr0_0.tsv.gz"
    df = pd.read_csv(path_file, sep='\t', compression='gzip')
    print(df.columns)
    df_select = df[['chr', col1, col2, 'counts_breast', 'counts_nonbreast',
                'p_value_X2_self', 'p_value_X2', 'p_value_F']]
    df_select = calculate_fc(df)
    df_select['info'] = df_select['chr'].map(str) + '_' + df_select[col1].map(str) + '_' + df_select[col2].map(str) + '_' + df_select['foldchange'].map(str) + '_' + df_select['bigger'].map(str)

    run_different_fc(df_select, type_file, non_coding, path_analyse, db, 'ALL')

    df_breast = df_select[df_select['bigger'] == 'b']
    run_different_fc(df_breast, type_file, non_coding, path_analyse, db, 'breast')

    df_nonbreast = df_select[df_select['bigger'] == 'nb']
    run_different_fc(df_nonbreast, type_file, non_coding, path_analyse, db, 'nonbreast')


def main():
    config = get_config('gearshift')
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/db_laatste_copy.db"
    db = Database(path_db)
    path_analyse = config['analyse_new'] #'D:/Hanze_Groningen/STAGE/test/' 

    # per_snp
    type_file = 'per_snp'
    for non_coding in ['ALL', 'Coding', 'NonCoding']: #'ALL', 'Coding', 'Non
        run_all_corrections(path_analyse, type_file, non_coding, db, 'pos_start', 'pos_end')

    # Region 1000 and 2000
    for type_file in ['Region_1000', 'Region_2000']: #'Region_1000', 
         for non_coding in ['ALL', 'Coding', 'NonCoding']: #'ALL', 'Coding', 'Non
             run_all_corrections(path_analyse, type_file, non_coding, db, 'start_region', 'stop_region')

    db.close()





if __name__ == '__main__':
    main()