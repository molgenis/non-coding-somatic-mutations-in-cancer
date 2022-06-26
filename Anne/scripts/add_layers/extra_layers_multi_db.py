#!/usr/bin/env python3

# Imports
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from collections import Counter


def close_to(db, chr, start_pos, end_pos, variant_file, donor_dict, donor_list, donor_cancer_list, filter_num,
             name_variant):
    """
    Selects all SNPs that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (variant_file).
    :param db:                The database object
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param variant_file:      The file after which the layers with their specific region are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :param filter_num:        The number for filtering the total read count. (the number must be equal to or
                              greater than)
    :param name_variant:      The name of the layer
    :return:
    """
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
                    WHERE snp.chr = '%s' AND snp.pos_start >= %s AND snp.pos_end <= %s AND 
                            sum_dosage_GT.total_read_count_sum > %s 
                            AND sum_dosage_GT.mutant_allele_read_count_sum > 0 AND snp.ID = sum_dosage_GT.snp_ID 
                            AND snp.%s = 1
                    GROUP BY sum_dosage_GT.snp_ID, sum_dosage_GT.donor_ID;
                    """ %
                      (str(chr), int(start_pos), int(end_pos), int(filter_num), str(name_variant)))

    results = db.cursor.fetchall()
    if len(results) > 0:
        # Loop over results
        for res in results:
            # Append donor_ID to donor_list_snp
            donor_list_snp.append(donor_dict[res['donor_ID']])
            # Append snp_ID to snp_ID_list
            snp_ID_list.append(res['snp_ID'])
            # Check if donor_ID already exist in the dictionary
            if donor_dict[res['donor_ID']] in donor_read_count:
                donor_read_count[donor_dict[res['donor_ID']]][0] += res['total_read_count_sum']
                donor_read_count[donor_dict[res['donor_ID']]][1] += res['dosages']
                donor_read_count[donor_dict[res['donor_ID']]][2] += 1
            else:
                donor_read_count[donor_dict[res['donor_ID']]] = [res['total_read_count_sum'], res['dosages'], 1]
        # Loop over dict: donor_read_count
        for key, value in donor_read_count.items():
            donor_index = donor_list.index(key)
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
        variant_file.write(str(filter_num) + '\t' + str(chr) + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_ID_list)) + '\t' + ','.join(map(str, snp_ID_list)) + '\t' + str(len(donor_list_snp)) + '\t' + str(
            donor_count) + '\t' + str(cancer_count) + '\n')
    else:
        # If no results (SNPs) are found in a certain region, a line will be written
        variant_file.write(str(filter_num) + '\t' + str(chr) + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_ID_list)) + '\t-\t-\t-\t-\n')


def make_file_extra_layers(df_variant, chr, name_variant, path_db, config, donor_dict, donor_list, donor_cancer_list):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (variant_file).
    :param df_variant:        Dataframe with only the chosen chromosome
    :param chr:               Chromosome number or letter (without chr)
    :param name_variant:      The name of the layer
    :param path_db:           The path to the database
    :param config:            Dictionary with as keys the name of the paths and as value the paths   
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return:
    """
    # Database connection
    db = Database(path_db)
    # The header for the file
    header_file = 'filter\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    # Make file
    f = open(f'{config["genes_eQTL_etc"]}{name_variant}_{chr}_num_snps_ALL_db.tsv', 'w')
    # Write file to file
    f.write(header_file)
    for index, row in df_variant.iterrows():
        close_to(db, row['#Chromosome'].replace('chr', ''), row['Start'], row['End'], f, donor_dict, donor_list,
                 donor_cancer_list, 33, name_variant)
    f.close()
    # Committing the current transactions
    db.mydb_connection.commit()


def main():
    # Call get_config
    config = get_config('gearshift')
    name_layer = sys.argv[1]
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    # Path to file with regions of the layer
    path_file = config[name_layer]
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    # Set arguments to list
    arg_multi_list = []
    # Loop over chromosomes
    for item in list(set(df_variant['#Chromosome'])):
        select_variant_df = df_variant.loc[df_variant['#Chromosome'] == item]
        arg_multi_list.append(
            (select_variant_df, item, name_layer, path_db, config, donor_dict, donor_list, donor_cancer_list))

    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=make_file_extra_layers, iterable=arg_multi_list)
    pool.close()
    pool.join()
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
