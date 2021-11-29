#!/usr/bin/env python3
import os.path
import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/35737116/runtimeerror-invalid-display-variable
plt.switch_backend('agg')
import sys

# Get path from input
path = sys.argv[1]
# Read the file
df = pd.read_csv(path, sep='\t', index_col=0)
# Get the basename of the file
basename = os.path.basename(path).split('.')[0]


def make_hist(data, col, basename, titles_dict):
    """

    :param data:
    :param col:
    :param basename:
    :param titles_dict:
    :return:
    """
    # Create a figure
    fig = plt.figure(figsize=(15, 15))
    # Create a histogram
    sns.histplot(data=data, x=col).set_title(f'{col}')
    # Give the plot a title (comes from the dictionary)
    fig.suptitle(f'{titles_dict[col]}', fontsize=20)
    # # Make the y-axis a logarithmic scale
    plt.yscale('log')
    plt.savefig(f'{sys.argv[2]}{basename}_{col}_seaborn_hist.png')
    # Save plot
    plt.clf()


def make_subplot(data, column, basename, titles_dict):
    """

    :param data:
    :param column:
    :param basename:
    :param titles_dict: 
    :return:
    """
    # Create a subplot figure, with the number of columns as vertical plots (one below the other).
    # So if the data frame contains 5 columns, 5 plots are plotted one below the other.
    figure, axes = plt.subplots(len(data.columns), 1, figsize=(20, 20))
    # How much there should be between the plots and rounds in terms of empty space
    figure.tight_layout(pad=8.0)
    # Walk across columns in the new dataframe (data)
    for i, col in enumerate(data.columns):
        # Create a histogram
        ax = sns.histplot(data=data, x=col, ax=axes[i])
        # Make the y-axis a logarithmic scale
        ax.set_yscale('log')
        # When the column name is AD, titles are created per plot
        if column == 'AD':
            # ax.set_xlabel('AD (Allelic depths)')
            if i == 0:
                ax.title.set_text("REF")
            else:
                ax.title.set_text(f"ALT{i}")
    # Give the entire plot a title (comes from the dictionary)
    figure.suptitle(f'{titles_dict[column]}', fontsize=20)
    # Save plot
    plt.savefig(
        f'{sys.argv[2]}{basename}_{column}_seaborn_hist.png')
    plt.clf()


# Dictionary with titles for the plots per column
titles_dict = {'AD': 'Allelic depths (AD)', 'AF': 'Allele fractions (AF)', 'DP': 'read depth (DP)',
               'F1R2': 'Count of reads in F1R2 pair orientation supporting each allele (F1R2)',
               'F2R1': 'Count of reads in F2R1 pair orientation supporting each allele (F2R1)', 'GT': 'Genotype (GT)',
               'PGT': 'Physical phasing haplotype information (PGT)',
               'PS': 'Phasing set (PS) (typically the position of the first variant in the set)',
               'SB': "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias. (SB)"}

# Loop across columns
for column in list(df.columns):
    print(column)
    # Remove all NA values from the dataframe in this column
    data = df[df[column].notna()]
    # Check if the column name is in the list.
    # When the column name appears in the list, it means that there are
    # multiple values in this column. These values are separated by commas
    if column in ['AD', 'AF', 'F2R1', 'F1R2', 'SB']:
        # Split the column by comma
        split_df = df[column].str.split(',', expand=True)
        # Create numeric values of all values from the new dataframe
        data = split_df.apply(pd.to_numeric)
        # Create subplots with histograms per column
        make_subplot(data, column, basename, titles_dict)
    elif column != 'PD':
        # If the column name is equal to DP, make all values numeric values
        if column in ['DP']:
            data[column] = pd.to_numeric(data[column], downcast="float")
        # Check if the data frame is not empty
        if not data.empty:
            # Create histogram of plot
            make_hist(data, column, basename, titles_dict)
        else:
            print(f'{column} EMPTY')
