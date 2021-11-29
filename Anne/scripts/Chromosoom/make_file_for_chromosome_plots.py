#!/usr/bin/env python3
import pandas as pd
import io
import sys
import os
# import seaborn as sns
# import matplotlib.pyplot as plt

# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    """

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

path = "D:/Hanze_Groningen/STAGE/VCF/dbSNP/"
path_file = "D:/Hanze_Groningen/STAGE/VCF/dbSNP/merge_manual_bowtie.vcf"
# path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/data/"
# path_file = f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/merge_vcf/dbSNP_filter/merge_manual_bowtie.vcf'
print(path_file)
# Get the basename of the file
basename = os.path.basename(path_file).split('.')[0]
# Read vcf file
df = read_vcf(path_file)#(sys.argv[1].strip())
df.reset_index(level=0, inplace=True)
df['index'] = 'SNP' + df['index'].astype(str)
print(df)


# def make_hist(data, col, basename, num, path):
#     """

#     :param data:
#     :param col:
#     :param basename:
#     :param titles_dict:
#     :return:
#     """
#     # Create a figure
#     fig = plt.figure(figsize=(15, 15))
#     # Create a histogram
#     sns.histplot(data=data, x=col).set_title(f'{col}')
#     # Give the plot a title (comes from the dictionary)
#     fig.suptitle(f'{num}', fontsize=20)
#     # # Make the y-axis a logarithmic scale
#     plt.yscale('log')
#     # /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/data/
#     # D:\Hanze_Groningen\STAGE\VCF\dbSNP\R\
#     plt.savefig(f'{path}{num}_{basename}.png')
#     # Save plot
#     plt.clf()


for tissue_num in list(df.columns)[10:]:
    new_df = df.iloc[:, 0:3]
    new2_df = new_df
    print(new2_df.columns)
    new_df['END_POS'] = new_df['POS']
    new_df[tissue_num] = df[tissue_num]
    remove_df = new_df[~new_df[tissue_num].str.contains(":.:.")]
    remove_df.drop(tissue_num, axis=1, inplace=True)
    print(tissue_num)
    if ':' in tissue_num:
        version = tissue_num.split(':')[0]
        num_tis = tissue_num.split('_')[1].split('.')[0]
        num = f'{num_tis}_v{version}'
    else:
        num = tissue_num.split('_')[1].split('.')[0]
    # print(num)
    # remove_df['Data'] = 100
    # print(remove_df)
    # print(len(remove_df))
    # make_hist(remove_df, 'POS', basename, num, path)
    # remove_df.to_csv(f'{path}{num}_{basename}.tsv', sep='\t', index=False, header=False)
    #'index', 'CHROM', 'POS'
    new2_df.rename(columns={'index': 'Sample', 'CHROM': 'Chr', 'POS': 'Start'}, inplace=True)
    new2_df.to_csv(f'{path}R_{num}_{basename}.tsv', sep='\t', index=False)
    


