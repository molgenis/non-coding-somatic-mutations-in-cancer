#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os
import glob

# Also takes the folder 1 higher, so that I can do the import after
sys.path.append("..")
from vcf_reader import read_vcf


# # https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
# def read_vcf(path):
#     """

#     :param path:
#     :return:
#     """
#     with open(path, 'r') as f:
#         lines = [l for l in f if not l.startswith('##')]
#     return pd.read_csv(
#         io.StringIO(''.join(lines)),
#         dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
#                'QUAL': str, 'FILTER': str, 'INFO': str},
#         sep='\t'
#     ).rename(columns={'#CHROM': 'CHROM'})


def make_plot_format_vcf(path, basename):
    basename = basename.split('.')[0]
    # Read vcf file
    df = read_vcf(path)#(sys.argv[1].strip())    
    print(df.head())
    print(df.columns)

def make_plot_format_other(path, basename, out_path):
    df = pd.read_csv(path, sep='\t')
    print(df.head())
    print(df.columns)
    select_columns = ['icgc_donor_id', 'project_code', 'icgc_sample_id', 'chromosome', 'chromosome_start',
                         'chromosome_end', 'assembly_version', 'reference_genome_allele', 'mutated_to_allele', 
                         'total_read_count', 'platform', 'sequencing_strategy'] ##'mutant_allele_read_count', 'gene_affected', 'transcript_affected',
    select_df = df[select_columns]
    select_df.rename(columns={'icgc_donor_id': 'donor_id', 'project_code': 'project_id', 'icgc_sample_id': 'tissue_id', 
                        'chromosome': 'CHROM', 'chromosome_start': 'FROM', 'chromosome_end': 'TO', 
                        'assembly_version': 'genome_version', 'reference_genome_allele': 'ref', 'mutated_to_allele': 'alt', 
                         'total_read_count': 'depth',  'platform': 'platform', 'sequencing_strategy': 'seq_strategy'}, inplace=True) 
   # print(select_df.head())

    # nan_values = select_df.isna()
    # nan_columns = nan_values.any()
    # columns_with_nan = select_df.columns[nan_columns].tolist()
    # print(columns_with_nan)
    select_df['depth'].fillna(0, inplace=True)
    select_df = select_df.astype({'FROM': 'int64', 'TO': 'int64', 'depth': 'int64'})
    select_df['CHROM'] = 'chr' + select_df['CHROM'].astype(str)
    select_df["ID"] = ""
    select_df.to_csv(f'{out_path}{basename}_db3.tsv', sep="\t", index=False)
    
        
# path = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/simple_somatic_mutation.open.BOCA-UK.tsv" 
# # path = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/merge_manual_bwa_aln.vcf"
path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/'
out_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/'
print(path)
path_files = f"{path}*.tsv"
for fname in glob.glob(path_files):
    print(f'fname {fname}')
    basename = os.path.basename(fname).split('%2F')[3]
    print(f'basename {basename}')
    print('------')
    type_file = 'xxx'
    if type_file == 'vcf':
        make_plot_format_vcf(fname, basename, out_path)
    else:
        make_plot_format_other(fname, basename, out_path)
# Get the basename of the file
# basename = os.path.basename(path) #.split('.')[0]
# type_file = 'xxx'
# if type_file == 'vcf':
#     make_plot_format_vcf(path, basename)
# else:
#     make_plot_format_other(path, basename)




