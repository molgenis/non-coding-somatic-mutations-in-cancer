#!/usr/bin/env python3
import os.path
import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt

def make_hist(data, col, basename, path, chrom):
    """

    :param data:
    :param col:
    :param basename:
    :param titles_dict:
    :return:
    """
    # Create a figure
    fig, ax = plt.subplots(figsize=(25, 15))
    # Create a histogram
    sns.histplot(data=data, x=col, bins=200).set_title(f'{col}')
    # Give the plot a title (comes from the dictionary)
    fig.suptitle(f'hoi', fontsize=20)
    # # Make the y-axis a logarithmic scale
    #plt.yscale('log')
    ax.set_xlim(chrom[0], chrom[1])
    plt.savefig(f'{path}0000{basename}_{col}_seaborn_hist.png')
    # Save plot
    plt.clf()
    print('plot')


print('test')
# Get path from input
path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/chromosoom/"#sys.argv[1]
# Read the file
df = pd.read_csv(f'{path}data/SS6004099_merge_manual_bowtie.tsv', sep='\t', names=['Sample', 'Chr', 'Start', 'END_POS', 'info'])
print(df.tail())
# Get the basename of the file
basename = os.path.basename(f'{path}data/SS6004099_merge_manual_bowtie.tsv').split('.')[0]

chrom = pd.read_csv(f'{path}length_chromosome.tsv', sep='\t', header=None, index_col=0)
chr22 = chrom.loc[['chr22']].values.tolist()[0]
make_hist(df, 'Start', basename, path, chr22)