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


def layer_run(df_variant, name_variant, df, path_save, chr):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file.
    :param df_variant:   The file with the regios of the layer
    :param name_variant: The name of the layer
    :param df:           Dataframe with only the chosen chromosome out of the database
    :param path_save:    Path to save the file
    :param chr:          Chromosome number or letter
    :return:
    """
    # The header for the file
    header_file = 'filter\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    # Make file
    f = open(f'{path_save}{name_variant}_chr{chr}_num_snps_ALL_file.tsv', 'w')
    f.write(header_file)
    # Loop over regions of the layer
    for index, row in df_variant.iterrows():
        # Filter df on chr, start and stop position
        filter_region = df[
            (df['chr'] == row['#Chromosome']) & (df['pos_start'] >= row['Start']) & (df['pos_end'] <= row['End'])]
        # Get list of donores
        donor_list = list(filter_region['donor_ID'])
        # Dictionary with donor_id as key and value as how often this donor_id occurs in donor_list
        donor_count = dict(Counter(donor_list))
        # Get list of cancers
        cancer_list = list(filter_region['cancer'])
        # Dictionary with cancer as key and value as how often this cancer occurs in cancer_list
        cancer_count = dict(Counter(cancer_list))
        if len(filter_region) > 0:
            # write to file
            f.write(str(33) + '\t' + str(row['#Chromosome']) + '\t' + str(row['Start']) + '\t' + str(
                row['End']) + '\t' + str(
                len(filter_region['snp_ID'])) + '\t' + ','.join(map(str, list(filter_region['snp_ID']))) + '\t' + str(
                len(donor_list)) + '\t' + str(
                donor_count) + '\t' + str(cancer_count) + '\n')
        else:
            # If no results (SNPs) are found in a certain region, a line will be written
            f.write(str(33) + '\t' + str(row['#Chromosome']) + '\t' + str(row['Start']) + '\t' + str(
                row['End']) + '\t-\t-\t-\t-\t-\n')

    f.close()


def run_code(config, df, chr, name_layer):
    """
    Reads file and filters file by chromosome. Calls the layer_run function.
    :param config:  Dictionary with as keys the name of the paths and as value the paths   
    :param df:      Dataframe with only the chosen chromosome out of the database
    :param chr:     Chromosome number or letter
    :param name_layer: The name of the layer
    :return:
    """
    # Get file with regions of the layer
    path_file = config[name_layer]
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    # Filter file on chromosome
    df_variant['#Chromosome'] = df_variant['#Chromosome'].str.replace('chr', '')
    path_save = config["genes_eQTL_etc"]
    # Call layer_run
    layer_run(df_variant, name_layer, df, path_save, chr)


def main():
    # Call get_config
    config = get_config('gearshift')
    name_layer = sys.argv[1]
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    # Make dataframe of information in the database
    df = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                            sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, 
                            snp.chr, snp.pos_start, snp.pos_end , snp.DNase, snp.TFBS, snp.UCNE
                    FROM project, sum_dosage_GT, snp 
                    WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                                sum_dosage_GT.donor_project_ID = project.ID AND 
                                (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                                sum_dosage_GT.total_read_count_sum >= 33 ;''', db.mydb_connection)
    chromosomes_list = list(set(df['chr']))

    cpus = mp.cpu_count()
    arg_multi_list = []
    # Loop over chromosomes
    for chrom in chromosomes_list:
        df_chr = df[df['chr'] == chrom]
        arg_multi_list.append((config, df_chr, chrom, name_layer))
    # Multiprocess
    pool = Pool(processes=cpus)
    pool.starmap(func=run_code, iterable=arg_multi_list)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
