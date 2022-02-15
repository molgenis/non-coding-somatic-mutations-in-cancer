#!/usr/bin/env python3
import pandas as pd
import sys
import os
import glob

# Also takes the folder 1 higher, so that I can do the import after
# sys.path.append("..")
sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf


# def make_plot_format_vcf(path, basename):
#     """
#
#     :param db:  the database object
#     :return:
#     """
#     basename = basename.split('.')[0]
#     # Read vcf file
#     df = read_vcf(path)  # (sys.argv[1].strip())
#     print(df.head())
#     print(df.columns)


def make_plot_format_other(path, basename, out_path):
    """
    Selected and renames columns and saved them
    :param path:        The path to the file
    :param basename:    The basename of the file
    :param out_path:    The path where the new data should be stored
    :return:
    """
    # Read the data
    df = pd.read_csv(path, sep='\t', compression='gzip')
    # List of selected columns
    select_columns = ['icgc_donor_id', 'project_code', 'icgc_specimen_id', 'icgc_sample_id', 'chromosome', 'chromosome_start',
                      'chromosome_end', 'assembly_version', 'reference_genome_allele', 'mutated_to_allele',
                      'total_read_count', 'platform',  'sequencing_strategy', 'mutant_allele_read_count']  ##'mutant_allele_read_count', 'gene_affected', 'transcript_affected',
    # Select the columns out of the dataframe
    select_df = df[select_columns]
    # Rename columns
    select_df.rename(columns={'icgc_donor_id': 'donor_id', 'project_code': 'project_id', 'icgc_specimen_id': 'specimen_id', 
                                'icgc_sample_id': 'icgc_sample_id',
                              'chromosome': 'chr', 'chromosome_start': 'pos_start', 'chromosome_end': 'pos_end',
                              'assembly_version': 'genome_version', 'reference_genome_allele': 'ref',
                              'mutated_to_allele': 'alt', 'total_read_count': 'total_read_count', 'platform': 'platform',
                              'sequencing_strategy': 'seq_strategy', 'mutant_allele_read_count':'mutant_allele_read_count'}, inplace=True)
    # If depth is not entered, fill in with 0. #TODO
    # select_df['depth'].fillna(0, inplace=True)
    # Change type or some columns
    select_df = select_df.astype({'pos_start': 'int64', 'pos_end': 'int64'}) #'depth': 'int64'
    # Add chr to the column CHROM. 1 > chr1 etc.
    # select_df['chr'] = 'chr' + select_df['CHROM'].astype(str)
    # Add column ID #TODO
    # select_df["ID"] = ""

    # replace field that's entirely space (or empty) with NULL
    # select_df = select_df['total_read_count'].replace(r'^\s*$', 'NULL', regex=True)
    # select_df = select_df['mutant_allele_read_count'].replace(r'^\s*$', 'NULL', regex=True)
    # Save dataframe
    select_df.to_csv(f'{out_path}{basename}_db_NEW.tsv.gz', sep="\t", index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})


def main():
    """

    """
    # The path to the data
    path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/'
    # The path where the new data should be stored
    out_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_data_db/'
    # download?fn=%2Fcurrent%2FProjects%2FALL-US%2Fsimple_somatic_mutation.open.ALL-US.tsv.gz
    path_files = f"{path}download*.tsv.gz"
    # Loop over all files in path that ends with .tsv
    for filename in glob.glob(path_files):
        # Basename of file
        basename = os.path.basename(filename).split('%2F')[3]
        print(basename)
        type_file = 'other'
        if type_file == 'vcf':
            print('hoi')
            # make_plot_format_vcf(filename, basename, out_path)
        else:
            # Call make_plot_format_other
            make_plot_format_other(filename, basename, out_path)


if __name__ == '__main__':
    main()
