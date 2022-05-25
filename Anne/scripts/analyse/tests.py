# import multiprocessing as mp
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats.distributions import chi2
from bioinfokit import analys, visuz
from scipy.stats import fisher_exact
import time
from scipy.special import factorial
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from fisher import pvalue_npy
from scipy.stats import chi2_contingency
from scipy.stats import uniform, randint


def chi_square_self(df, colname_b, colname_nb, num_donor_b, num_donor_nb):
    """
    
    """
    df['counts_breast_plus'] = df[colname_b] + (1/num_donor_b)
    print(set(df['counts_breast_plus']))
    df['counts_nonbreast_plus'] = df[colname_nb] + (1/num_donor_nb)
    print(set(df['counts_nonbreast_plus']))
    n = np.array([num_donor_b, num_donor_nb])
    print(n)
    S = df[['counts_breast_plus', 'counts_nonbreast_plus']].to_numpy()
    f = n / n.sum()
    E1 = S.sum(axis=1)[:,None]*f
    C = (n-S)
    print(set(C[0]))
    print(set(C[1]))
    E2 = C.sum(axis=1)[:,None]*f
    D1 = ((E1-S)**2/E1).sum(axis=1)
    D2 = ((E2-C)/E2).sum(axis=1)
    X2 = D1 + D2
    #
    p_value_X2 = chi2.sf(X2,1)
    #
    log10_p_value = -np.log10(p_value_X2)
    log10_p_value
    #
    df['X2'] = X2
    df['p_value_X2_self'] = p_value_X2
    df['log10_p_value_X2_self'] = log10_p_value
    #
    S_C = np.concatenate((S, C), axis=1)
    return df, S_C, n, S, C    

def chi_square(S_C, df):
    """
    
    """
    start_time = time.perf_counter()

    p_value_X2_TEST = list()
    for index, value in enumerate(S_C):
        table = np.array([[value[0], value[2]], [value[1], value[3]]])
        stat, p, dof, expected = chi2_contingency(table)
        p_value_X2_TEST.append(p)
        if (index%1000000) == 0:
            print(len(p_value_X2_TEST))
    print(len(p_value_X2_TEST))

    print("--- %s seconds ---" % (time.perf_counter() - start_time))
    df['p_value_X2'] = p_value_X2_TEST
    return df
    

def log2_fc(df, n, S):
    """
    
    """
    constant = np.log(n[0]) - np.log(n[1])
    log2_fc = (np.log(S[:, 0]) - np.log(S[:, 1]) - constant) / np.log(2)
    df['log2_fc'] = log2_fc # log2(FC)
    return df



def fisher_test(S_C, df):
    """
    
    """
    start_time = time.perf_counter()

    p_value_F = list()
    for index, value in enumerate(S_C):
        table = np.array([[value[0], value[2]], [value[1], value[3]]])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        p_value_F.append(p)
        if (index%1000000) == 0:
            print(len(p_value_F))
    print(len(p_value_F))

    print("--- %s seconds ---" % (time.perf_counter() - start_time))
    df['p_value_F'] = p_value_F
    return df


def shapiro_test(df, column_name, column_name_filter, path_file, type_analyse, type_df):
    """
    As the p value obtained from the Shapiro-Wilk test is significant (p < 0.05), 
    we conclude that the data is not normally distributed. Further, in histogram data 
    distribution shape does not look normal. Therefore, Mann-Whitney U test is more appropriate for analyzing two samples.
    """
    w, pvalue = stats.shapiro(df[column_name])
    print(f'w:{w}, pvalue:{pvalue}')
    # plot histogram
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Frequency histogram')
    ax1.hist(df[column_name_filter], bins=df[column_name].max(), histtype='bar', ec='k') 
    ax2.hist(df[column_name], bins=df[column_name].max(), histtype='bar', ec='k') 
    plt.title(f"{type_analyse} - {type_df} nonbreast")
    plt.savefig(f'{path_file}shapiro_test_histo_{type_analyse}_{type_df}_{column_name}.png')
    plt.clf()


def mannwhitney(df, column_breast, column_nonbreast):
    """
    Mann-Whitney U test interpretation: As the p value obtained from the Mann-Whitney U test 
    is significant (U = 489.5, p < 0.05), we conclude that the yield of the two genotypes 
    significantly different from each other .
    """
    # U1, p = mannwhitneyu(list(breast['filter_snps']), list(nonbreast['filter_snps']), alternative = 'two-sided')
    print(mannwhitneyu(list(df[column_breast]), list(df[column_nonbreast]), alternative = 'two-sided'))



def volcano_plot(df, p_value_column, path_file, type_analyse, type_df):
    """
    
    """
    visuz.GeneExpression.volcano(df=df,lfc='log2_fc',pv=p_value_column, plotlegend=True, legendpos='upper right')
    plt.title(f"{type_analyse} - {type_df} nonbreast")
    plt.savefig(f'{path_file}volcano_{p_value_column}_{type_analyse}_{type_df}_both.png')
    plt.clf()
    

def all_test(df, num_donor_b, num_donor_nb, type_df, type_analyse, path_file):
    """
    
    """
    print('\nboxplot')
    df.boxplot(column=['counts_breast'], grid=False)
    plt.title(f"{type_analyse} - {type_df} Breast")
    plt.savefig(f'{path_file}boxplot_{type_analyse}_{type_df}_breast.png')
    plt.clf()
    df.boxplot(column=['counts_nonbreast'], grid=False)
    plt.title(f"{type_analyse} - {type_df} nonbreast")
    plt.savefig(f'{path_file}boxplot_{type_analyse}_{type_df}_nonbreast.png')
    plt.clf()
    
    print('\nfilter columns (divide by max)')
    df['filter_snps_b'] = df['counts_breast']/ num_donor_b
    df['filter_snps_nb'] = df['counts_nonbreast']/ num_donor_nb
    
    # print('\nshapiro_test')
    # shapiro_test(df, 'counts_breast', 'filter_snps_b', path_file, type_analyse, type_df)
    # shapiro_test(df, 'counts_nonbreast', 'filter_snps_nb', path_file, type_analyse, type_df)
    
    # print('\nmannwhitney')
    # mannwhitney(df, 'filter_snps_b', 'filter_snps_nb')
    
    print('\ntests')
    df, S_C, n, S, C = chi_square_self(df, 'counts_breast', 'counts_nonbreast', num_donor_b, num_donor_nb)
    df = chi_square(S_C, df)
    df = log2_fc(df, n, S)
    df = fisher_test(S_C, df)
    df = log2_fc(df, n, S)
    df.to_csv(f"{path_file}{type_analyse}_{type_df}_both_0_TESTS.tsv", sep='\t', encoding='utf-8', index=False)
    
    # print('\nvolcano_plot')
    # volcano_plot(df, 'p_value_X2_self', path_file, type_analyse, type_df)
    # volcano_plot(df, 'p_value_X2', path_file, type_analyse, type_df)
    # volcano_plot(df, 'p_value_F', path_file, type_analyse, type_df)
    
    return df


