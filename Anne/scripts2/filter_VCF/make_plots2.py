#!/usr/bin/env python3
import os.path
import pandas as pd
import io
import csv
import seaborn as sns
import matplotlib.pyplot as plt
# https://stackoverflow.com/questions/35737116/runtimeerror-invalid-display-variable
plt.switch_backend('agg')
import sys

def format_plots(column, col, sample_name, df, dir_figure):
    """
    
    :param column: 
    :param sample_name: 
    :param df: 
    :param dir_figure: 
    :return: 
    """
    data = df[df[column].notna()]    
    if col in ['AD', 'AF', 'F2R1', 'F1R2', 'SB','MBQ', 'MMQ', 'NALOD', 'NLOD', 'POPAF', 'TLOD']:
        # Split the column by comma
        split_df = data[column].str.split(',', expand=True)
        # Create numeric values of all values from the new dataframe
        data = split_df.apply(pd.to_numeric)
        # Create subplots with histograms per column
        make_subplot(data, column, col, sample_name, dir_figure)  # , titles_dict)
    elif col  not in ['PS', 'MFRL', 'MPOS', 'RPA']:
        # If the column name is equal to DP, make all values numeric values
        if column in ['DP',    'GERMQ', 'PON', 'SEQQ', 'STR', 'STRANDQ', 'STRQ']:
            data[column] = pd.to_numeric(data[column], downcast="float")
        # Check if the data frame is not empty
        if not data.empty:
            # Create histogram of plot
            make_hist(data, column, col, sample_name, dir_figure)  # , titles_dict)
        else:
            print(f'{column} EMPTY')


def make_hist(data, column, col, sample_name, dir_figure):  # , titles_dict):
    """

    :param data:
    :param col:
    :param sample_name:
    :param dir_figure:
    :return:
    """
    # Create a figure
    fig = plt.figure(figsize=(15, 15))
    # Create a histogram
    sns.histplot(data=data, x=column).set_title(f'{col}')
    # Give the plot a title (comes from the dictionary)
    fig.suptitle(f'{sample_name}-{col}', fontsize=20)
    # # Make the y-axis a logarithmic scale
    plt.yscale('log')
    plt.savefig(f'{dir_figure}{sample_name}_{col}_seaborn_hist.png')
    # Save plot
    plt.clf()
    plt.close()


def make_subplot(data, column, col, sample_name, dir_figure):
    """

    :param data:
    :param column:
    :param sample_name:
    :param dir_figure: 
    :return:
    """
    # Create a subplot figure, with the number of columns as vertical plots (one below the other).
    # So if the data frame contains 5 columns, 5 plots are plotted one below the other.
    figure, axes = plt.subplots(len(data.columns), 1, figsize=(20, 20))
    # How much there should be between the plots and rounds in terms of empty space
    figure.tight_layout(pad=8.0)
    # Walk across columns in the new dataframe (data)
    for i, col_df in enumerate(data.columns):
        # Create a histogram
        ax = sns.histplot(data=data, x=col_df, ax=axes[i])
        # Make the y-axis a logarithmic scale
        ax.set_yscale('log')
        # When the column name is AD, titles are created per plot
        if col == 'AD':
            # ax.set_xlabel('AD (Allelic depths)')
            if i == 0:
                ax.title.set_text("REF")
            else:
                ax.title.set_text(f"ALT{i}")
    # Give the entire plot a title (comes from the dictionary)
    figure.suptitle(f'{sample_name}-{col}', fontsize=20)
    # Save plot
    plt.savefig(
        f'{dir_figure}{sample_name}_{col}_seaborn_hist.png')
    plt.clf()
    plt.close()



# Get path from input
# path = "D:/Hanze_Groningen/STAGE/VCF/aln_SS6004099.tsv"#sys.argv[1]
# path = "D:/Hanze_Groningen/STAGE/VCF/mem_SS6005042-mem_SS6005043.vcf"
path = sys.argv[1]
print(path)
# Read the file
# https://stackoverflow.com/questions/18016037/pandas-parsererror-eof-character-when-reading-multiple-csv-files-to-hdf5
df = pd.read_csv(path, sep='\t', quoting=csv.QUOTE_NONE, encoding='utf-8')
# Get the basename of the file
# basename = os.path.basename(path).split('.')[0]
df_info_format = df.iloc[:, list(df.columns).index("[7]FILTER"):]

for column in list(df_info_format.columns):
    print(column)
    path_to, filename = os.path.split(path)
    filename, extension = os.path.splitext(filename)
    dir_figure=sys.argv[2]
    if filename.split('-')[0] in column:
        sample_name = filename.split('-')[0]
        col = column.split(':')[-1]
        format_plots(column, col, sample_name, df_info_format, dir_figure)
    elif len(filename.split('-'))>1 and filename.split('-')[1] in column:
        sample_name = filename.split('-')[1]
        col = column.split(':')[-1]
        format_plots(column, col, sample_name, df_info_format, dir_figure)
    else:
        col = column.split(']')[1]
        format_plots(column, col, filename, df_info_format, dir_figure)
