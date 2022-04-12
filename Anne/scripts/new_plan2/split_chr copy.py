from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix
from multiprocessing import Pool, Queue
import multiprocessing as mp
import math
from search_snps_between import close_to, write_sparse_matrix



def loop_over_genes(db, chr, length_chrom, steps, donor_dict, donor_list, donor_cancer_list,
                        save_path, region_list, sparse_matrix_region, sparse_matrix_overall, filter_num, with_type):
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
    # TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    # TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    # TODO
    # The header for the files before_gene_file and after_gene_file
    header_file = 'filter_num\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    # Make file
    before_gene_file = open(f'{save_path}chr{chr}_region_{steps}.tsv', 'w')
    # Write header
    before_gene_file.write(header_file)
    # TODO make nist
    total_read = [0] * len(donor_list)
    # Loop over genes in file
    i_start = 0
    for i in range(0, length_chrom+1, steps):
        if i != 0:
            print(i_start, '-', i)
            # Call close_to
            sparse_matrix_region, sparse_matrix_overall, region_list, total_read = close_to(db,
                                                                                            f'{chr}:{i_start}-{i}',
                                                                                            chr,
                                                                                            i_start,
                                                                                            i,
                                                                                            before_gene_file,
                                                                                            donor_dict,
                                                                                            donor_list,
                                                                                            region_list,
                                                                                            sparse_matrix_region,
                                                                                            sparse_matrix_overall,
                                                                                            donor_cancer_list,
                                                                                            total_read, filter_num, with_type)
        if length_chrom == (i + (length_chrom % steps)):
            last_i = (i + (length_chrom % steps))
            if i != last_i:
                print(i+1, '-', last_i)
                sparse_matrix_region, sparse_matrix_overall, region_list, total_read = close_to(db,
                                                                                                f'{chr}:{i_start}-{i}',
                                                                                                chr,
                                                                                                (i + 1),
                                                                                                last_i,
                                                                                                before_gene_file,
                                                                                                donor_dict,
                                                                                                donor_list,
                                                                                                region_list,
                                                                                                sparse_matrix_region,
                                                                                                sparse_matrix_overall,
                                                                                                donor_cancer_list,
                                                                                                total_read, filter_num, with_type)
        i_start = (i + 1)
    # Close file
    before_gene_file.close()
    # Call write_sparse_matrix
    write_sparse_matrix(sparse_matrix_region, region_list, donor_list, save_path, 'region',
                        donor_cancer_list, total_read, chr)    
    write_sparse_matrix(sparse_matrix_overall, region_list, donor_list, save_path, 'overall',
                        donor_cancer_list, total_read, chr)

def set_region_list(steps, chr, length_chrom):
    """
    
    :param db:  The database object
    :return:
    """
    # name_file = f'D:/Hanze_Groningen/STAGE/NEW PLAN/chr{key}.tsv.gz'
    i_start=0
    region_list = []
    for i in range(0,length_chrom+1,steps):
        if i != 0:
            region_list.append(f'{chr}:{i_start}-{i}')
            if length_chrom == (i + (length_chrom % steps)):
                last_i = (i + (length_chrom % steps))
                if i != last_i:
                    region_list.append(f'{chr}:{(i + 1)}-{last_i}')
        i_start = (i + 1)
    return region_list


def multiprocess(chr_length, steps, donor_list, donor_dict, donor_cancer_list, save_path, filter_num, with_type):
    """
    
    :param db:  The database object
    :return:
    """
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"  # /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    for chr, length_chrom in chr_length.items():
        region_list = set_region_list(steps, chr, length_chrom)
        # Creating a len(donor_list) * len(gene_name_list) sparse matrix
        sparse_matrix = csr_matrix((len(donor_list), len(region_list)),
                                                dtype=np.int8).toarray()
        print('loop over genes')
        # Call loop_over_genes
        loop_over_genes(db, chr, length_chrom, steps, donor_dict, donor_list, donor_cancer_list,
                        save_path, region_list, sparse_matrix, sparse_matrix.copy(), filter_num, with_type)
        print('CLOSE')

def split_dict(d, n, chr_length):
    """
    # https://stackoverflow.com/questions/22878743/how-to-split-dictionary-into-multiple-dictionaries-fast

    
    :param db:  The database object
    :return:
    """
    keys = list(d.keys())
    for i in range(0, len(keys), n):
        yield {k: chr_length[k] for k in keys[i: i + n]}


def main():
    # Path of the database
    path_db = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydb_L.db"  # /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    # Path to save files
    save_path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/UMAP/"
    # The steps for the region
    steps= 2000

    #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-247
    #https://www.nature.com/articles/nature07517#Sec6
    # SNP discovery increases with increasing depth: essentially all homozygous positions are detected at 15×, whereas heterozygous positions accumulate more gradually to 33× (Fig. 5a). 
    filter_num = 33
    with_type = 'chrom'
    #https://en.wikipedia.org/wiki/Human_genome
    chr_length = {'1':248956422, '2':242193529, '3':198295559, '4':190214555,
                '5':181538259, '6':170805979, '7':159345973, '8':145138636,
                '9':138394717, '10':133797422, '11':135086622, '12':133275309,
                '13':114364328, '14':107043718, '15':101991189, '16':90338345,
                '17':83257441, '18':80373285, '19':58617616, '20':64444167,
                '21':46709983, '22':50818468, 'X':156040895, 'Y':57227415} #{'17':83257441, '18':80373285, '19':58617616, '20':64444167}
    
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)


    cpus = mp.cpu_count()-2
    n = math.ceil(len(chr_length) / cpus)
    print(n)
    print(len(chr_length))

    arg_multi_list = []
    for item in split_dict({i: i for i in chr_length}, n, chr_length):
        print(item)
        arg_multi_list.append((item, steps, donor_list, donor_dict, donor_cancer_list, save_path, filter_num, with_type))

    pool = Pool(processes=cpus)
    pool.starmap(func=multiprocess, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # for chr, length_chrom in chr_length.items():
    #     region_list = set_region_list(steps, chr, length_chrom)
    #     print(region_list)
    #     # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    #     sparse_matrix = csr_matrix((len(donor_list), len(region_list)),
    #                                             dtype=np.int8).toarray()
    #     print('loop over genes')
    #     # Call loop_over_genes
    #     loop_over_genes(db, chr, length_chrom, steps, donor_dict, donor_list, donor_cancer_list,
    #                     save_path, region_list, sparse_matrix, sparse_matrix.copy())
    #     print('CLOSE')
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
