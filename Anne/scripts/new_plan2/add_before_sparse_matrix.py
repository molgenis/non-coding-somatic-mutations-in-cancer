#!/usr/bin/env python3
from Database import Database
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import sys

sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
import glob

def merge_sparce_matrix(save_path, whole_numpy_array):
    """

    :param save_path: string with path of where the file should be saved.
    :param whole_numpy_array: empty sparce array of len(donor_list) at len(gene_name_list)+1
    :return:
    """
    path = f"{save_path}*_sparsematrix_bef_overall.tsv.gz"
    for fname in glob.glob(path):
        print(fname)
        df = pd.read_csv(fname, sep='\t', compression='gzip')
        cancer_list = df['cancer']
        donor_id = df['donor_id']
        columns_name = df.columns
        df.set_index('donor_id', inplace=True)
        df.drop(['cancer'], axis=1, inplace=True)
        numpy_array = df.to_numpy()
        whole_numpy_array = whole_numpy_array + numpy_array
    # Make dataframe of numpy array
    df_whole = pd.DataFrame(data=whole_numpy_array, columns=columns_name)
    # Add donor IDs as column
    df_whole['donor_id'] = donor_id
    # Add cancer types as column
    df_whole['cancer'] = cancer_list
    # Makes the donor_ids column the index of the data frame
    df_whole.set_index('donor_id', inplace=True)
    df_whole.to_csv(f'{save_path}ALL_sparsematrix_bef_overall.tsv.gz', sep="\t", index=True, encoding='utf-8',
              compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})

def main():
    config = get_config()
    # Path of the database
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    # Gene file path
    gene_path = config['all_genes'] 
    # Path to save files
    save_path = config['umap_path'] 
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    # Gene name list
    gene_name_list = gene_df['hg19.kgXref.geneSymbol'].tolist()
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    # +1 for the column total_reads
    whole_numpy_array = csr_matrix((len(donor_list), len(gene_name_list)+1),
                                             dtype=np.int8).toarray() 
    # Call merge_sparce_matrix
    merge_sparce_matrix(save_path, whole_numpy_array)

if __name__ == '__main__':
    main()