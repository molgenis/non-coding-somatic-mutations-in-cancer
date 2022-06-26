#!/usr/bin/env python3

#Imports

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def make_plot(new_df, cancer_type, type_plot, path_file):
    """
    Make boxplots
    :param new_df: Dataframe
    :param cancer_type: Type of cancer
    :param type_plot: Type of plot
    :param path_file: Path in which the files are saved
    :return:
    """
    new_df.plot(kind='box')
    plt.title("Number of snps in region of 1000000bp")
    plt.savefig(f'{path_file}box_ALL_1.png')
    plt.clf()

    new_df[cancer_type].plot(kind=type_plot)
    plt.title("Number of snps in region of 1000000bp")
    plt.savefig(f'{path_file}{type_plot}_{cancer_type}_1.png')
    plt.clf()

    ax = sns.boxplot(data=new_df).set(ylabel='#snps')
    plt.savefig(f'{path_file}box_ALL_2.png')
    plt.clf()

    if type_plot == 'box':
        ax = sns.boxplot(data=new_df[cancer_type], palette="Set2").set(xlabel=cancer_type, ylabel='#snps',
                                                                       title='Number of snps in region of 1000000bp')
    else:
        ax = sns.displot(new_df[cancer_type]).set(xlabel=cancer_type, title='Number of snps in region of 1000000bp')
    plt.savefig(f'{path_file}{type_plot}_{cancer_type}_2.png')
    plt.clf()


def outliers(df, column, cancer_type, path_file):
    """
    Get outliers
    :param df: Dataframe
    :param column: Name of the column
    :param cancer_type: Type of cancer
    :param path_file: Path in which the files are saved
    :results: outliers: dataframe with the outliers
    """
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    print(
        f'Q1: {Q1}, Q3: {Q3}, IQR: {IQR}, #Outlier: {((df[column] < (Q1 - 1.5 * IQR)) | (df[column] > (Q3 + 1.5 * IQR))).sum()}')
    outliers = df[((df[column] < (Q1 - 1.5 * IQR)) | (df[column] > (Q3 + 1.5 * IQR)))]
    outliers.sort_values(by=[1, 2], ascending=[False, True], inplace=True)
    print(outliers)
    header = ['num', 'snps', 'chr', 'start_region', 'end_region']
    outliers.to_csv(f'{path_file}outliers_{cancer_type}.tsv', sep='\t', encoding='utf-8', index=False, header=header)
    outliers['INFO'] = outliers[2].astype(str) + '_' + outliers[3].astype(str) + ':' + outliers[4].astype(str)

    return outliers


def main():
    path_file = 'D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/'
    # nonbreast cancer
    nonbreast = pd.read_csv(f'{path_file}nonbreast_LONG_ALL.tsv', sep='\t', header=None)
    outliers_nb = outliers(nonbreast, 1, 'nonbreast', path_file)
    # breast cancer
    breast = pd.read_csv(f'{path_file}breast_LONG_ALL.tsv', sep='\t', header=None)
    outliers_b = outliers(breast, 1, 'breast', path_file)
    # See overlap and unique outliers
    match = list(set(outliers_nb['INFO']) & set(outliers_b['INFO']))
    unique_nb = list(set(outliers_nb['INFO']) - set(outliers_b['INFO']))
    unique_b = list(set(outliers_b['INFO']) - set(outliers_nb['INFO']))
    # Write the unqiue outliers for breast cancer to a file
    unique_breast = outliers_b[outliers_b['INFO'].isin(unique_b)]
    header = ['num', 'snps', 'chr', 'start_region', 'end_region']
    unique_breast.drop('INFO', axis=1, inplace=True)
    unique_breast.to_csv(f'{path_file}outliers_unique_breast.tsv', sep='\t', encoding='utf-8', index=False,
                         header=header)
    print(f'Tot_breast: {len(outliers_b)}, Tot_nonbreast: {len(outliers_nb)}')
    print(f'match: {len(match)}, unique_breast: {len(unique_b)}, unique_nonbreast: {len(unique_nb)}')


if __name__ == '__main__':
    main()
