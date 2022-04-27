# Import libraries
from turtle import title
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import diff
import seaborn as sns
import pandas as pd

def make_plot(new_df, cancer_type, type_plot, path_file):
    new_df.plot(kind='box')
    # plt.ylabel('#snps')
    plt.title("Number of snps in region of 1000000bp")
    plt.savefig(f'{path_file}box_ALL_1.png')
    plt.clf()


    new_df[cancer_type].plot(kind=type_plot)
    # plt.ylabel('#snps')
    plt.title("Number of snps in region of 1000000bp")
    plt.savefig(f'{path_file}{type_plot}_{cancer_type}_1.png')
    plt.clf()


    # sns.set_theme(style="whitegrid")
    ax = sns.boxplot(data=new_df).set(ylabel='#snps')
    plt.savefig(f'{path_file}box_ALL_2.png')
    plt.clf()

    if type_plot == 'box':
        ax = sns.boxplot(data=new_df[cancer_type], palette="Set2").set(xlabel=cancer_type, ylabel='#snps', title='Number of snps in region of 1000000bp')
    else:
        ax = sns.displot(new_df[cancer_type]).set(xlabel=cancer_type, title='Number of snps in region of 1000000bp')
    # plt.show()
    plt.savefig(f'{path_file}{type_plot}_{cancer_type}_2.png')
    plt.clf()


def outliers(df, column, cancer_type, path_file):
    
    print(cancer_type)
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    print(f'Q1: {Q1}, Q3: {Q3}, IQR: {IQR}, #Outlier: {((df[column] < (Q1 - 1.5 * IQR)) | (df[column] > (Q3 + 1.5 * IQR))).sum()}')
    outliers = df[((df[column]<(Q1-1.5*IQR)) | (df[column]>(Q3+1.5*IQR)))]
    outliers.sort_values(by=[1, 2], ascending=[False, True], inplace=True)
    print(outliers)
    header = ['num', 'snps', 'chr', 'start_region', 'end_region']
    outliers.to_csv(f'{path_file}outliers_{cancer_type}.tsv', sep='\t', encoding='utf-8', index=False, header=header)
    outliers['INFO'] = outliers[2].astype(str) + '_' + outliers[3].astype(str) + ':' + outliers[4].astype(str) 

    return outliers



def main():
    path_file = 'D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/' 
    #"D:/Hanze_Groningen/STAGE/R/PLOTS/kary/vs/before/chr3/breast_LONG_chr3.tsv"

    nonbreast = pd.read_csv(f'{path_file}nonbreast_LONG_ALL.tsv', sep='\t', header=None)
    outliers_nb = outliers(nonbreast, 1, 'nonbreast', path_file)

    breast = pd.read_csv(f'{path_file}breast_LONG_ALL.tsv', sep='\t', header=None)
    outliers_b = outliers(breast, 1, 'breast', path_file)

    match = list(set(outliers_nb['INFO']) & set(outliers_b['INFO']))
    unique_nb = list(set(outliers_nb['INFO']) - set(outliers_b['INFO']))
    unique_b = list(set(outliers_b['INFO']) - set(outliers_nb['INFO']))

    unique_breast = outliers_b[outliers_b['INFO'].isin(unique_b)]
    header = ['num', 'snps', 'chr', 'start_region', 'end_region']
    unique_breast.drop('INFO', axis=1, inplace=True)
    unique_breast.to_csv(f'{path_file}outliers_unique_breast.tsv', sep='\t', encoding='utf-8', index=False, header=header)
    print(f'Tot_breast: {len(outliers_b)}, Tot_nonbreast: {len(outliers_nb)}')
    print(f'match: {len(match)}, unique_breast: {len(unique_b)}, unique_nonbreast: {len(unique_nb)}')




    # new_df = nonbreast[[1]]
    # new_df.rename(columns={ new_df.columns[0]: "nonbreast" }, inplace = True)
    # new_df['breast'] = breast[1]
    # # BOX
    # cancer_type = 'breast' # breast nonbreast
    # type_plot = 'box' # hist box
    # make_plot(new_df, cancer_type, type_plot, path_file)    
    # cancer_type = 'nonbreast' # breast nonbreast
    # make_plot(new_df, cancer_type, type_plot, path_file)
    # # HIST
    # cancer_type = 'breast' # breast nonbreast
    # type_plot = 'hist' # hist box
    # make_plot(new_df, cancer_type, type_plot, path_file)
    # cancer_type = 'nonbreast' # breast nonbreast
    # make_plot(new_df, cancer_type, type_plot, path_file)

    

if __name__ == '__main__':
    main()
