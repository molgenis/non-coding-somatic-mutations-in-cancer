#!/usr/bin/env python3

#Imports
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.distributions import chi2
from bioinfokit import analys, visuz
from scipy.stats import fisher_exact
import time
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency
from scipy.stats.contingency import relative_risk


def chi_square_self(df, colname_b, colname_nb, num_b, num_nb):
    """
    Chi-square test performed manually
    :param df: The dataframe with the breast and nonbreast counts
    :param colname_b:  Column name for breast cancer counts
    :param colname_nb: Column name for nonbreast cancer counts
    :param num_b:  Number of breast cancer donors or SNPs
    :param num_nb:  Number of nonbreast cancer donors or SNPs
    :return:    df: Dataframe with new columns
                S_C: array with counts_breast_plus and counts_nonbreast_plus and Total number of donors/SNPs min
                     breast cancer or non breast cancer donors/SNPs
                n: array with num_b and num_donor_nb
                S: array with counts_breast_plus and counts_nonbreast_plus
                C:  Total number of donors/SNPs min breast cancer or non breast cancer donors/SNPs
    """
    # All plus 1/n (where n is the number of donors or snps)
    df['counts_breast_plus'] = df[colname_b] + (1 / num_b)
    df['counts_nonbreast_plus'] = df[colname_nb] + (1 / num_nb)
    # Total number of donors
    n = np.array([num_b, num_nb])
    S = df[['counts_breast_plus', 'counts_nonbreast_plus']].to_numpy()
    f = n / n.sum()
    E1 = S.sum(axis=1)[:, None] * f
    # Total number of donors/SNPs min breast cancer or non breast cancer donors/SNPs
    C = (n - S)
    E2 = C.sum(axis=1)[:, None] * f
    D1 = ((E1 - S) ** 2 / E1).sum(axis=1)
    D2 = ((E2 - C) / E2).sum(axis=1)
    X2 = D1 + D2
    # p_value
    p_value_X2 = chi2.sf(X2, 1)
    # log10 p_value
    log10_p_value = -np.log10(p_value_X2)
    log10_p_value
    # Add new columns to df
    df['X2'] = X2
    df['p_value_X2_self'] = p_value_X2
    df['log10_p_value_X2_self'] = log10_p_value
    # Merge S and C
    S_C = np.concatenate((S, C), axis=1)
    return df, S_C, n, S, C


def chi_square(S_C, df):
    """
    Running chi-square
    :param S_C: array with counts_breast_plus and counts_nonbreast_plus and Total number of donors/SNPs min
                breast cancer or non breast cancer donors/SNPs
    :param df:    The dataframe with the breast and nonbreast counts   
    :return:  df:  Dataframe with new columns
    """
    # Keep track of how long it takes
    start_time = time.perf_counter()
    # List of p_values
    p_value_X2 = list()
    # Loop over S_C
    for index, value in enumerate(S_C):
        table = np.array([[value[0], value[2]], [value[1], value[3]]])
        stat, p, dof, expected = chi2_contingency(table)
        # Append p-value to list
        p_value_X2.append(p)

    print("--- %s seconds ---" % (time.perf_counter() - start_time))
    df['p_value_X2'] = p_value_X2
    return df


def log2_fc(df, n, S):
    """
    Calculate log2 fold change
    :param df: The dataframe with the breast and nonbreast counts  
    :param n:  array with num_b and num_donor_nb
    :param S: array with counts_breast_plus and counts_nonbreast_plus
    :return:  df:  Dataframe with new columns
    """
    constant = np.log(n[0]) - np.log(n[1])
    log2_fc = (np.log(S[:, 0]) - np.log(S[:, 1]) - constant) / np.log(2)
    df['log2_fc'] = log2_fc
    return df


def fisher_test(S_C, df):
    """
    Running the fisher exact test
    :param S_C: array with counts_breast_plus and counts_nonbreast_plus and Total number of donors/SNPs min
                breast cancer or non breast cancer donors/SNPs
    :param df:   The dataframe with the breast and nonbreast counts      
    :return: df:    Dataframe with new columns
    """
    # Keep track of how long it takes
    start_time = time.perf_counter()
    # List with p-values
    p_value_F = list()
    # Loop over S_C
    for index, value in enumerate(S_C):
        table = np.array([[value[0], value[2]], [value[1], value[3]]])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        # Append p-value to list
        p_value_F.append(p)

    print("--- %s seconds ---" % (time.perf_counter() - start_time))
    df['p_value_F'] = p_value_F
    return df


def shapiro_test(df, column_name, column_name_filter, path_file, type_analyse, type_df):
    """
    Running the Shapiro test
    As the p value obtained from the Shapiro-Wilk test is significant (p < 0.05), 
    we conclude that the data is not normally distributed. Further, in histogram data 
    distribution shape does not look normal. Therefore, Mann-Whitney U test is more appropriate for analyzing two samples.
    :param df: The dataframe with the breast and nonbreast counts
    :param column_name:  column name
    :param column_name_filter: column name of the filter columns
    :param path_file: The path to which the image will be written       
    :param type_analyse:   name of the layer 
    :param type_df:        noncoding, coding or mix 
    :return:    
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
    Running the MannWhitney test
    Mann-Whitney U test interpretation: As the p value obtained from the Mann-Whitney U test 
    is significant (U = 489.5, p < 0.05), we conclude that the yield of the two genotypes 
    significantly different from each other .
    :param df: The dataframe with the breast and nonbreast counts
    :param column_breast:  column name for breast cancer counts
    :param column_nonbreast: column name for nonbreast cancer counts
    :return:    
    """
    print(mannwhitneyu(list(df[column_breast]), list(df[column_nonbreast]), alternative='two-sided'))


def volcano_plot(df, p_value_column, path_file, type_analyse, type_df):
    """
    Make volcano plots
    :param df: The dataframe with the breast and nonbreast counts
    :param p_value_column:  The name of the p-value column
    :param path_file: The path to which the image will be written
    :param type_analyse:  name of the layer 
    :param type_df:        noncoding, coding or mix 
    :return:    
    """
    visuz.GeneExpression.volcano(df=df, lfc='log2_fc', pv=p_value_column, plotlegend=True, legendpos='upper right')
    plt.title(f"{type_analyse} - {type_df} nonbreast")
    plt.savefig(f'{path_file}volcano_{p_value_column}_{type_analyse}_{type_df}_both.png')
    plt.clf()


def calculate_relative_risk(S_C, df):
    """
    Running relative risk test
    :param S_C: array with counts_breast_plus and counts_nonbreast_plus and Total number of donors/SNPs min breast
                cancer or non breast cancer donors/SNPs
    :param df:      The dataframe with the breast and nonbreast counts 
    :return:   df: Dataframe with new columns
    """
    relative_risk_values = list()
    confidence_interval = list()
    confidence_interval_bon = list()
    one_in_interval = list()
    one_in_interval_bon = list()
    list_values_RR = list()
    # Loop over S_C
    for index, value in enumerate(S_C):
        # Run relative risk
        result = relative_risk(int(round(value[0])), int(round(value[0] + value[2])), int(round(value[1])),
                               int(round(value[1] + value[3])))
        list_values_RR.append([value[0], (value[0] + value[2]), value[1], (value[1] + value[3])])
        relative_risk_values.append(result.relative_risk)
        # Get confidence interval
        con_interval = result.confidence_interval(confidence_level=0.95)
        confidence_interval.append(f'{con_interval[0]}-{con_interval[1]}')
        # See if 1 is within the confidence interval
        if con_interval[0] < 1 and con_interval[1] > 1:
            # If 1 is not in it, it is not significant
            one_in_interval.append(False)
        else:
            # When 1 is in it it is significant
            one_in_interval.append(True)
        # Confidence interval with bonferroni correction
        con_interval_bon = result.confidence_interval(confidence_level=1 - (0.05 / len(df)))
        confidence_interval_bon.append(f'{con_interval_bon[0]}-{con_interval_bon[1]}')
        if con_interval_bon[0] < 1 and con_interval_bon[1] > 1:
            # If 1 is not in it, it is not significant
            one_in_interval_bon.append(False)
        else:
            # When 1 is in it it is significant
            one_in_interval_bon.append(True)

    # Add new columsn to dataframe
    df['relative_risk_values'] = relative_risk_values
    df['confidence_interval'] = confidence_interval
    df['one_in_interval'] = one_in_interval
    df['confidence_interval_bon'] = confidence_interval_bon
    df['one_in_interval_bon'] = one_in_interval_bon
    df['list_values_RR'] = list_values_RR
    return df


def all_test(df, num_b, num_nb, type_df, type_analyse, path_file, select_chrom, i):
    """
    Calls all tests
    :param df: The dataframe with the breast and nonbreast counts 
    :param num_donor_b:  Number of breast cancer donors or SNPs
    :param num_nb:  Number of nonbreast cancer donors or SNPs
    :param type_df:        noncoding, coding or mix 
    :param type_analyse:  name of the layer 
    :param path_file: The path to which the image will be written
    :param select_chrom:  When multiprocessing is used it indicates which chromosome is selected.
                          When no multiprocessing is used it is always chr0.
    :param i:       When multiprocessing is used, the i indicates in which part it is currently working.
                    When no multiprocessing is used, i is always 0.
    :return:    df: Dataframe with new columns
    """
    print('\nboxplot')
    df.boxplot(column=['counts_breast'], grid=False)
    plt.title(f"{type_analyse} - {type_df} Breast")
    plt.savefig(f'{path_file}boxplot_{type_analyse}_{type_df}_breast_{select_chrom}_{i}.png')
    plt.clf()
    df.boxplot(column=['counts_nonbreast'], grid=False)
    plt.title(f"{type_analyse} - {type_df} nonbreast")
    plt.savefig(f'{path_file}boxplot_{type_analyse}_{type_df}_nonbreast_{select_chrom}_{i}.png')
    plt.clf()

    print('\nfilter columns (divide by max)')
    df['filter_snps_b'] = df['counts_breast'] / num_b
    df['filter_snps_nb'] = df['counts_nonbreast'] / num_nb

    print('\nshapiro_test')
    shapiro_test(df, 'counts_breast', 'filter_snps_b', path_file, type_analyse, type_df)
    shapiro_test(df, 'counts_nonbreast', 'filter_snps_nb', path_file, type_analyse, type_df)

    print('\nmannwhitney')
    mannwhitney(df, 'filter_snps_b', 'filter_snps_nb')

    print('\ntests')
    df, S_C, n, S, C = chi_square_self(df, 'counts_breast', 'counts_nonbreast', num_b, num_nb)
    df = chi_square(S_C, df)
    df = log2_fc(df, n, S)
    df = fisher_test(S_C, df)
    df = log2_fc(df, n, S)
    df.to_csv(f"{path_file}new/{type_analyse}_{type_df}_both_0_TESTS_{select_chrom}_{i}_NO_RR.tsv", sep='\t',
              encoding='utf-8', index=False)
    df = calculate_relative_risk(S_C, df)
    df.to_csv(f"{path_file}new/{type_analyse}_{type_df}_both_0_TESTS_{select_chrom}_{i}.tsv", sep='\t',
              encoding='utf-8', index=False)

    print('\nvolcano_plot')
    volcano_plot(df, 'p_value_X2_self', path_file, type_analyse, type_df)
    volcano_plot(df, 'p_value_X2', path_file, type_analyse, type_df)
    volcano_plot(df, 'p_value_F', path_file, type_analyse, type_df)

    return df
