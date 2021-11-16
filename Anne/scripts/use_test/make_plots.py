#!/usr/bin/env python3
# importing os.path module
import os.path
import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/35737116/runtimeerror-invalid-display-variable
plt.switch_backend('agg')
import sys

path = sys.argv[1]
df = pd.read_csv(path, sep='\t', index_col=0)
basename = os.path.basename(path).split('.')[0]
print(basename)
print(df.head())


def make_hist(data, column, basename, titles_dict):
    fig = plt.figure(figsize=(15, 15))
    sns.histplot(data=data, x=column).set_title(f'{column}')
    fig.suptitle(f'{titles_dict[column]}', fontsize=20)
    plt.yscale('log')
    plt.savefig(f'{sys.argv[2]}{basename}_{column}_seaborn_hist.png')
    plt.clf()


titles_dict = {'AD': 'Allelic depths (AD)', 'AF': 'Allele fractions (AF)', 'DP': 'read depth (DP)',
              'F1R2': 'Count of reads in F1R2 pair orientation supporting each allele (F1R2)',
              'F2R1': 'Count of reads in F2R1 pair orientation supporting each allele (F2R1)', 'GT': 'Genotype (GT)',
              'PGT': 'Physical phasing haplotype information (PGT)',
              'PS': 'Phasing set (PS) (typically the position of the first variant in the set)',
              'SB': "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias. (SB)"}

for column in list(df.columns):
    print(column)
    data = df[df[column].notna()]
    if column in ['AD', 'AF', 'F2R1', 'F1R2', 'SB']:
        split_df = df[column].str.split(',', expand=True)
        data = split_df.apply(pd.to_numeric)
        f, axes = plt.subplots(len(data.columns), 1, figsize=(20, 20))
        f.tight_layout(pad=8.0)
        for i, col in enumerate(data.columns):
            if not data.empty:
                ax = sns.histplot(data=data, x=col, ax=axes[i])
                ax.set_yscale('log')
                if column == 'AD':
                    ax.set_xlabel('AD (Allelic depths)')
                    if i == 0:
                        ax.title.set_text("REF")
                    else:
                        ax.title.set_text(f"ALT{i}")
            else:
                print(f'{col} EMPTY')
        f.suptitle(f'{titles_dict[column]}', fontsize=20)
        plt.savefig(
            f'{sys.argv[2]}{basename}_{column}_seaborn_hist.png')
        plt.clf()
    elif column == 'PID':
        print(column)
    else:
        if column in ['DP']:
            data[column] = pd.to_numeric(data[column], downcast="float")
        if not data.empty:
            make_hist(data, column, basename, titles_dict)
        else:
            print(f'{column} EMPTY')
