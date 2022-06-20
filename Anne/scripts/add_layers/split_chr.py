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
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



def check_snp_id(db, snp_ids, donor_dict, donor_list, gene_index, sparse_matrix_overall, sparse_matrix_region,
                 total_read):
    """
    Checks for each snp which donors have this snp. 
    For this snp, these donors must meet total_read_count > 0 and mutant_allele_read_count > 0.
    These donors are then added to a list. 
    This list will eventually contain the donors (a donor can be more frequent) who contain snps in a certain region.
    :param db:                The database object
    :param snp_ids:           A list of snp IDs that fall within a particular region.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :return: donor_list_snp : List of donors (research donor IDs) who contain snps within a particular region.
    """
    #
    donor_read_count = dict()
    # Make donor_list_snp. List because double donors count
    donor_list_snp = list()
    # Loop over list with snp IDs
    for snp_id in snp_ids:
        # Make donor_set. Set because duplicate snps from a donor don't count
        donor_set = set()
        # See how many snp_IDs there are in this table that are equal to the ID of the snp and total_read_count > 0
        # AND mutant_allele_read_count >
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID, total_read_count, dosages
                        FROM 'donor_has_snp'
                        WHERE snp_ID = %s AND total_read_count > 0 AND mutant_allele_read_count > 0;
                        """ % snp_id)
        results = db.cursor.fetchall()
        for res in results:
            # Add donor_ID to the set
            donor_set.add(donor_dict[res['donor_ID']])  # res['donor_ID']
            if donor_dict[res['donor_ID']] in donor_read_count:
                donor_read_count[donor_dict[res['donor_ID']]][0] += res['total_read_count']
                donor_read_count[donor_dict[res['donor_ID']]][1] += res['dosages']
                donor_read_count[donor_dict[res['donor_ID']]][2] += 1
            else:
                donor_read_count[donor_dict[res['donor_ID']]] = [res['total_read_count'], res['dosages'], 1]
        # Extend the list with the donor_set
        donor_list_snp.extend(donor_set)
    for key, value in donor_read_count.items():
        donor_index = donor_list.index(key)
        sparse_matrix_overall[donor_index, gene_index] = value[2]
        # donor_read_count[key] = value[2] / (value[1] / value[0])
        sparse_matrix_region[donor_index, gene_index] = value[2] / (value[1] / value[0])
        if (value[2] / (value[1] / value[0])) <= 0:
            print(f'{value[2]} / ({value[1]} / {value[0]})')
        total_read[donor_index] += value[0]

    return donor_list_snp, sparse_matrix_overall, sparse_matrix_region, total_read


def close_to(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, region_list, sparse_matrix_region,
             sparse_matrix_overall, donor_cancer_list, total_read):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (gene_file).
    :param db:                The database object
    :param gene:              The name of the gene
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param gene_file:         The file after which the genes with their specific region before or after
                              a gene are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param region_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param sparse_matrix_region:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return: sparse_matrix_region:   A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
             region_list:  List of gene names (to be used later as columns in the sparse matrix)
    """
    # Get index of the gene name in the list
    gene_index = region_list.index(gene)
    snp_id_list = db.get_snps(chr, start_pos, end_pos)
    if len(snp_id_list) > 0:
        # Call donor_list_before
        donor_list_snp, sparse_matrix_overall, sparse_matrix_region, total_read = check_snp_id(db, snp_id_list, donor_dict,
                                                                                            donor_list, gene_index,
                                                                                            sparse_matrix_overall,
                                                                                            sparse_matrix_region,
                                                                                            total_read)
        # Make cancer_list
        cancer_list = list()
        for donor in donor_list_snp:
            donor_index = donor_list.index(donor)
            # Adds +1 to the sparse_matrix at the position of donor region
            # sparse_matrix_overall[donor_index, gene_index] += 1
            cancer_list.append(donor_cancer_list[donor_index])

        # Creates a dictionary from the list with as key the name (cancer type or donor ID) and
        # as value how often that name occurs in the list.
        cancer_count = dict(Counter(cancer_list))
        donor_count = dict(Counter(donor_list_snp))

        # Write to file
        gene_file.write(chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_id_list)) + '\t' + ','.join(map(str, snp_id_list)) + '\t' + str(len(donor_list_snp)) + '\t' + str(
            donor_count) + '\t' + str(cancer_count) + '\n')
        # str(len(set(donor_list_snp)))+'\t'+','.join(map(str,donor_list_snp))+'\t'
    else:
        # Write to file
        gene_file.write(chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_id_list)) + '\t-\t-\t-\t-\n')

    return sparse_matrix_region, sparse_matrix_overall, region_list, total_read


def write_sparse_matrix(sparse_matrix, region_list, donor_list, save_path, pos, donor_cancer_list, total_read, chr):
    """
    Writes the sparse_matrix to a compressed (.gz) .tsv file
    :param sparse_matrix:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param region_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param save_path:         Path to save files
    :param pos:               A string indicating whether it is before or after a gene (bef/aft)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return:

    """
    # Creates a dataframe from the sparse_matrix with column names region_list
    df = pd.DataFrame(data=sparse_matrix, columns=region_list)
    # Add donor IDs as column
    df['donor_id'] = donor_list
    # Add cancer types as column
    df['cancer'] = donor_cancer_list
    # Add total read count as column
    df['total_read'] = total_read
    # Makes the donor_ids column the index of the data frame
    df.set_index('donor_id', inplace=True)
    # Write the dataframe to a compressed .tsv file
    df.to_csv(f'{save_path}chr{chr}_sparsematrix_{pos}.tsv.gz', sep="\t", index=True, encoding='utf-8',
              compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})


def loop_over_genes(db, chr, length_chrom, steps, donor_dict, donor_list, donor_cancer_list,
                        save_path, region_list, sparse_matrix_region, sparse_matrix_overall):
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
    header_file = 'chr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
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
                                                                                            total_read)
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
                                                                                                total_read)
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


def multiprocess(chr_length, steps, donor_list, donor_dict, donor_cancer_list, save_path, config):
    """
    
    :param db:  The database object
    :return:
    """
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"  
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
                        save_path, region_list, sparse_matrix, sparse_matrix.copy())
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
    config = get_config('gearshift')
    # Path of the database
    path_db = config['database'] 
    # Database connection
    db = Database(path_db)
    # Path to save files
    save_path = config['umap_path'] 
    # The steps for the region
    steps= 2000
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
        arg_multi_list.append((item, steps, donor_list, donor_dict, donor_cancer_list, save_path, config))

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
