from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix
# from matplotlib import pyplot as plt
from search_snps_between import close_to, write_sparse_matrix
import sys
from multiprocessing import Pool, Queue
import multiprocessing as mp

sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
import glob


def main():
    config = get_config()
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" 
    # Database connection
    db = Database(path_db)
    
    # check_filter(db)
    # Path of the genes and there positions #'/local/1161112/rawdata/cancer_data/genes_eQTL_etc/all_genes_new.tsv'

    gene_path = config['all_genes'] # snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
    # Path to save files
    save_path = config['umap_path'] #TODO umap_path: '/local/1161112/rawdata/cancer_data/UMAP/'

    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    # print(len(gene_df))
    # # Replace all empty values with NaN in the column proteinID
    # gene_df['proteinID'].replace('', np.nan, inplace=True)
    # # Drop all NaN values (in column proteinID)
    # gene_df.dropna(subset=['proteinID'], inplace=True)
    # print(len(gene_df))
    gene_name_list = gene_df['hg19.kgXref.geneSymbol'].tolist()
    print(len(gene_name_list))
    """
    Okosun et al. 
    * Regions of  -2000bp - 250bp (5' UTRs if applicable) from the transcription starting sites (TSS) 
    for each transcript were screened. For transcripts from the same gene that share the same promoter 
    mutation profiles, only one representative transcript was selected.
    """
    # Call add_value
    # add_value(db)
    print('set GENE')
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    whole_numpy_array = csr_matrix((len(donor_list), len(gene_name_list)+1),
                                             dtype=np.int8).toarray() #+1 for total_reads
    # df_whole = pd.DataFrame(index=range(len(donor_list)), columns = gene_name_list)
    # df_whole['donor_list'] = donor_list
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
        # whole_numpy_array = np.add(whole_numpy_array, numpy_array) 
        # df_whole.add(df)
    df_whole = pd.DataFrame(data=whole_numpy_array, columns=columns_name)
    # Add donor IDs as column
    df_whole['donor_id'] = donor_id
    # Add cancer types as column
    df_whole['cancer'] = cancer_list
    # Makes the donor_ids column the index of the data frame
    df_whole.set_index('donor_id', inplace=True)
    df_whole.to_csv(f'{save_path}ALL_sparsematrix_bef_overall.tsv.gz', sep="\t", index=True, encoding='utf-8',
              compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})


    # # Creates a dataframe from the sparse_matrix with column names gene_name_list
    # df = pd.DataFrame(data=sparse_matrix, columns=gene_name_list)
    # # Add donor IDs as column
    # df['donor_id'] = donor_list
    # # Add cancer types as column
    # df['cancer'] = donor_cancer_list
    # # Add total read count as column
    # df['total_read'] = total_read
    # # Makes the donor_ids column the index of the data frame
    # df.set_index('donor_id', inplace=True)
    # # Write the dataframe to a compressed .tsv file
    # df.to_csv(f'{save_path}{part_num}_sparsematrix_{pos}.tsv.gz', sep="\t", index=True, encoding='utf-8',
    #           compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})



if __name__ == '__main__':
    main()