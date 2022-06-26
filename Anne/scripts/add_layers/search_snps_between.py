#!/usr/bin/env python3

# Imports
import pandas as pd
from collections import Counter


def close_to(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, gene_name_list,
             sparse_matrix_overall, donor_cancer_list, total_read, filter_num, with_type):
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
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param sparse_matrix_overall: A matrix which contains very few non-zero elements. It contains the counts of a
                                specific region-donor combination.
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :param total_read:          Lift of total read counts
    :param filter_num:          The number for filtering the total read count. (the number must be equal to or
                                greater than)
    :param with_type:           Check if it is gene or not
    :return: sparse_matrix_overall: A matrix which contains very few non-zero elements. It contains the counts of a
                                   specific region-donor combination.
             gene_name_list:  List of gene names (to be used later as columns in the sparse matrix)
             total_read:          Lift of total read counts             
    """
    # Get index of the gene name in the list
    gene_index = gene_name_list.index(gene)
    if with_type == 'genes':
        # Replace the name with gene_name:chromosoom:start_pos-end_pos
        gene_name_list[gene_index] = f'{gene}:chr{chr}:{int(start_pos)}-{int(end_pos)}'
    # Dictionary with the key donor_id and as value a list with the total number of read counts, 
    # total dosage and how often this donor occurs
    donor_read_count = dict()
    # Make donor_list_snp. List because double donors count
    donor_list_snp = list()
    # List of SNPs
    snp_ID_list = list()

    db.cursor.execute("""
                    SELECT sum_dosage_GT.snp_ID , sum_dosage_GT.donor_ID, 
                           sum_dosage_GT.total_read_count_sum , sum_dosage_GT.mutant_allele_read_count_sum,
                           sum_dosage_GT.dosages
                    FROM snp, sum_dosage_GT
                    WHERE snp.chr = '%s' AND snp.pos_start >= %s AND snp.pos_end <= %s 
                    AND sum_dosage_GT.total_read_count_sum > %s AND sum_dosage_GT.mutant_allele_read_count_sum > 0 
                    AND (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND snp.ID = sum_dosage_GT.snp_ID
                    GROUP BY sum_dosage_GT.snp_ID, sum_dosage_GT.donor_ID;
                    """ %
                      (str(chr), int(start_pos), int(end_pos), int(filter_num)))

    results = db.cursor.fetchall()
    if len(results) > 0:
        # Loop over results
        for res in results:
            # Append donor_ID to donor_list_snp
            donor_list_snp.append(donor_dict[res['donor_ID']])
            # Append snp_ID to snp_ID_list
            snp_ID_list.append(res['snp_ID'])
            # Check if donor_ID already exists
            if donor_dict[res['donor_ID']] in donor_read_count:
                donor_read_count[donor_dict[res['donor_ID']]][0] += res['total_read_count_sum']
                donor_read_count[donor_dict[res['donor_ID']]][1] += res['dosages']
                donor_read_count[donor_dict[res['donor_ID']]][2] += 1
            else:
                donor_read_count[donor_dict[res['donor_ID']]] = [res['total_read_count_sum'], res['dosages'], 1]
        # Loop over dict: donor_read_count
        for key, value in donor_read_count.items():
            donor_index = donor_list.index(key)
            sparse_matrix_overall[donor_index, gene_index] = value[2]
            total_read[donor_index] += value[0]
        # Make cancer_list
        cancer_list = list()
        for donor in donor_list_snp:
            donor_index = donor_list.index(donor)
            cancer_list.append(donor_cancer_list[donor_index])
        # Creates a dictionary from the list with as key the name (cancer type or donor ID) and
        # as value how often that name occurs in the list.
        cancer_count = dict(Counter(cancer_list))
        donor_count = dict(Counter(donor_list_snp))
        # Write to file
        if with_type == 'genes':
            gene_file.write(
                str(filter_num) + '\t' + gene + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
                    len(snp_ID_list)) + '\t' + ','.join(map(str, snp_ID_list)) + '\t' + str(
                    len(donor_list_snp)) + '\t' + str(
                    donor_count) + '\t' + str(cancer_count) + '\n')
        else:
            gene_file.write(filter_num + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
                len(snp_ID_list)) + '\t' + ','.join(map(str, snp_ID_list)) + '\t' + str(
                len(donor_list_snp)) + '\t' + str(
                donor_count) + '\t' + str(cancer_count) + '\n')
    else:
        # If no results (SNPs) are found in a certain region, a line will be written
        if with_type == 'genes':
            gene_file.write(
                str(filter_num) + '\t' + gene + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
                    len(snp_ID_list)) + '\t-\t-\t-\t-\n')
        else:
            gene_file.write(filter_num + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
                len(snp_ID_list)) + '\t-\t-\t-\t-\n')
    return sparse_matrix_overall, gene_name_list, total_read


def write_sparse_matrix(sparse_matrix, gene_name_list, donor_list, save_path, pos, donor_cancer_list, total_read,
                        part_num):
    """
    Writes the sparse_matrix to a compressed (.gz) .tsv file
    :param sparse_matrix:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param save_path:         Path to save files
    :param pos:               A string indicating whether it is before or after a gene (bef/aft)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :param total_read:          Lift of total read counts
    :param part_num:            The number of the part for multiprocessing
    :return:

    """
    # Creates a dataframe from the sparse_matrix with column names gene_name_list
    df = pd.DataFrame(data=sparse_matrix, columns=gene_name_list)
    # Add donor IDs as column
    df['donor_id'] = donor_list
    # Add cancer types as column
    df['cancer'] = donor_cancer_list
    # Add total read count as column
    df['total_read'] = total_read
    # Makes the donor_ids column the index of the data frame
    df.set_index('donor_id', inplace=True)
    # Write the dataframe to a compressed .tsv file
    df.to_csv(f'{save_path}{part_num}_sparsematrix_{pos}_NEW.tsv.gz', sep="\t", index=True, encoding='utf-8',
              compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
