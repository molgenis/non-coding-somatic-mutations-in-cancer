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



def loop_over_genes(config, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list,
                    save_path, gene_name_bef, gene_name_aft, sparse_matrix_before_region, sparse_matrix_after_region,
                    sparse_matrix_before_overall, sparse_matrix_after_overall, filter_num, with_type, part_num):
    """
    Loop over the gene data frame.
    :param db:                  The database object
    :param gene_df:             The data frame with all genes (columns: #name, chrom, strand, txStart, txEnd, cdsStart,
                                cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param position_out_gene:   Region before the start position of a gene or after the stop position of a gene
    :param position_in_gene:    Region after the start position of a gene or before the stop position of a gene
    :param donor_dict:          A dictionary with as key the automatically generated donor ID and as value the donor
                                IDs that are used in the research.
    :param donor_list:          List of donor names (to be used later as rows in the sparse matrix)
    :param donor_cancer_list:   List of cancers. This list has the same order as donor_list.
    :param save_path:           Path to save files
    :param gene_name_bef:       List of donor names before (to be used later as rows in the sparse matrix)
    :param gene_name_aft:       List of donor names after (to be used later as rows in the sparse matrix)
    :param sparse_matrix_before_region:A matrix which contains very few non-zero elements. It contains the counts of a
                                specific region (before gene)-donor combination.
    :param sparse_matrix_after_region: A matrix which contains very few non-zero elements. It contains the counts of a
                                specific region (after gene)-donor combination.
    :param sparse_matrix_before_overall:
    :param sparse_matrix_after_overall:
    :return:

    """
    print(part_num)
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" 
    # Database connection
    db = Database(path_db)
     # The header for the files before_gene_file and after_gene_file
    header_file = 'filter\tgene\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    # Make file
    before_gene_file = open(f'{save_path}{part_num}_gene_before_{position_out_gene}_{position_in_gene}.tsv', 'w')
    after_gene_file = open(f'{save_path}{part_num}_gene_after_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    before_gene_file.write(header_file)
    after_gene_file.write(header_file)
    # TODO make nist
    total_read_before = [0] * len(donor_list)
    total_read_after = [0] * len(donor_list)
    # Loop over genes in file
    for index, row in gene_df.iterrows():
        # Remove 'chr' from the chromosome (chr1 --> 1)
        chr = row['hg19.knownGene.chrom'].replace('chr', '')
        # print(chr, row['hg19.knownGene.txStart'], row['hg19.knownGene.txEnd'])
        # Call close_to
        print('BEF')
        sparse_matrix_before_region, sparse_matrix_before_overall, gene_name_bef, total_read_before = close_to(db,
                                                                                                        row['hg19.kgXref.geneSymbol'],
                                                                                                        chr,
                                                                                                        row[
                                                                                                            'hg19.knownGene.txStart'] - position_out_gene,
                                                                                                        row[
                                                                                                            'hg19.knownGene.txStart'] + position_in_gene,
                                                                                                        before_gene_file,
                                                                                                        donor_dict,
                                                                                                        donor_list,
                                                                                                        gene_name_bef,
                                                                                                        sparse_matrix_before_region,
                                                                                                        sparse_matrix_before_overall,
                                                                                                        donor_cancer_list,
                                                                                                        total_read_before, filter_num, with_type)
        print('AFT')
        sparse_matrix_after_region, sparse_matrix_after_overall, gene_name_aft, total_read_after = close_to(db, row['hg19.kgXref.geneSymbol'],
                                                                                                      chr,
                                                                                                      row[
                                                                                                          'hg19.knownGene.txEnd'] - position_in_gene,
                                                                                                      row[
                                                                                                          'hg19.knownGene.txEnd'] + position_out_gene,
                                                                                                      after_gene_file,
                                                                                                      donor_dict,
                                                                                                      donor_list,
                                                                                                      gene_name_aft,
                                                                                                      sparse_matrix_after_region,
                                                                                                      sparse_matrix_after_overall,
                                                                                                      donor_cancer_list,
                                                                                                      total_read_after, filter_num, with_type)

    # Close file
    before_gene_file.close()
    after_gene_file.close()
    # Call write_sparse_matrix
    write_sparse_matrix(sparse_matrix_before_region, gene_name_bef, donor_list, save_path, 'bef_region',
                        donor_cancer_list, total_read_before, part_num)
    write_sparse_matrix(sparse_matrix_after_region, gene_name_aft, donor_list, save_path, 'aft_region',
                        donor_cancer_list, total_read_after, part_num)
    write_sparse_matrix(sparse_matrix_before_overall, gene_name_bef, donor_list, save_path, 'bef_overall',
                        donor_cancer_list, total_read_before, part_num)
    write_sparse_matrix(sparse_matrix_after_overall, gene_name_aft, donor_list, save_path, 'aft_overall',
                        donor_cancer_list, total_read_after, part_num)


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
    # Region before the start position of a gene or after the stop position of a gene
    position_out_gene = 2000
    # Region after the start position of a gene or before the stop position of a gene
    position_in_gene = 250
    #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-247
    #https://www.nature.com/articles/nature07517#Sec6
    # SNP discovery increases with increasing depth: essentially all homozygous positions are detected at 15×, whereas heterozygous positions accumulate more gradually to 33× (Fig. 5a). 
    filter_num = 33
    with_type = 'genes'
    # Call add_value
    # add_value(db)
    print('set GENE')
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    sparse_matrix_before_region = csr_matrix((len(donor_list), len(gene_name_list)),
                                             dtype=np.int8).toarray()
    sparse_matrix_after_region = csr_matrix((len(donor_list), len(gene_name_list)),
                                            dtype=np.int8).toarray()



    cpus = mp.cpu_count()

    list_part_num = list()
    arg_multi_list = []
    for i in range(0, len(gene_df)+101, 100):
        df = gene_df.iloc[i:i+100]
        list_part_num.append(i)
        arg_multi_list.append((config, df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list,
                    save_path, gene_name_list, gene_name_list.copy(), sparse_matrix_before_region,
                    sparse_matrix_after_region, sparse_matrix_before_region.copy(), sparse_matrix_after_region.copy(), filter_num, with_type, i))

    print(len(list_part_num))
    print(list_part_num)

    pool = Pool(processes=cpus)
    pool.starmap(func=loop_over_genes, iterable=arg_multi_list)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()

