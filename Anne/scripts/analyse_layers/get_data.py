#!/usr/bin/env python3

#Imports
import sys
import pandas as pd

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def get_data_db(filter_par, path_file, path_db):
    """
    Retrieves information from the database
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore
                       uses the previously stored data.
    :param path_file: Path where the files are saved
    :param path_db: Path to the database
    :return: df_whole: Dataframe with all the information out of the database
    """
    if filter_par:
        # Database connection
        db = Database(path_db)
        df_whole = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                                    sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, sum_dosage_GT.GT2,
                                    snp.chr, snp.pos_start, snp.pos_end, snp.ref, snp.alt,
                                    snp.UCNE, snp.TFBS, snp.DNase, 
                                    snp.in_transcript, snp.in_coding, snp.in_exon,
                                    snp.before_gene, snp.after_gene, 
                                    snp.ID_eQTL, snp.eQTL, snp.close_eQTL_3000
                            FROM project, sum_dosage_GT, snp 
                            WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                                    sum_dosage_GT.donor_project_ID = project.ID AND 
                                    (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                                    sum_dosage_GT.total_read_count_sum >= 33;''', db.mydb_connection)
        db.close()
        df_whole.to_csv(f"{path_file}df_whole.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        df_whole = pd.read_csv(f"{path_file}df_whole.tsv", sep='\t')

    df_whole = df_whole.loc[df_whole['chr'] != 'MT']
    return df_whole


def filter_whole_cancer(filter_par, df_whole, path_file):
    """
    Filter or grab the breast cancer data and the nonbreast cancer data
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param df_whole: Dataframe with all the information out of the database
    :param path_file:  Path where the files are saved
    :return: all_breast: All the SNPs of the breast cancer donors
             all_nonbreast: All the SNPs of the nonbreast cancer donors
    """
    if filter_par:
        # Filters out breast cancer data and not breast cancer data
        all_breast = df_whole.loc[df_whole['cancer'] == 'Breast']
        all_breast.to_csv(f"{path_file}all_breast.tsv", sep='\t', encoding='utf-8', index=False)
        all_nonbreast = df_whole.loc[df_whole['cancer'] != 'Breast']
        all_nonbreast.to_csv(f"{path_file}all_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        # Get the data frames for breast cancer and not breast cancer
        all_breast = pd.read_csv(f"{path_file}all_breast.tsv", sep='\t')
        all_nonbreast = pd.read_csv(f"{path_file}all_nonbreast.tsv", sep='\t')
    return all_breast, all_nonbreast


def filter_noncoding(filter_par, df_whole, path_file):
    """
    Filter or grab the noncoding data
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param df_whole: Dataframe with all the information out of the database
    :param path_file:  Path where the files are saved
    :return: noncoding_df: Dataframe with the noncoding information (SNPs)
    """
    if filter_par:
        # Get the noncoding SNPs
        noncoding_df = df_whole.loc[
            (df_whole['in_transcript'] == 0) & (df_whole['in_coding'] == 0) & (df_whole['in_exon'] == 0)]
        noncoding_df.to_csv(f"{path_file}noncoding_df.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        # Get the dataframe for the noncoding SNPs
        noncoding_df = pd.read_csv(f"{path_file}noncoding_df.tsv", sep='\t')
    return noncoding_df


def filter_noncoding_cancer(filter_par, noncoding_df, path_file):
    """
    Filter the noncoding data on breast cancer and nonbreast cancer
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param noncoding_df: Dataframe with the noncoding information (SNPs)
    :param path_file:  Path where the files are saved
    :return: noncoding_breast: All the noncoding SNPs of the breast cancer donors
             noncoding_nonbreast All the noncoding SNPs of the nonbreast cancer donors
    """
    if filter_par:
        # Filters out breast cancer data and not breast cancer data
        noncoding_breast = noncoding_df.loc[noncoding_df['cancer'] == 'Breast']
        noncoding_breast.to_csv(f"{path_file}noncoding_breast.tsv", sep='\t', encoding='utf-8', index=False)
        noncoding_nonbreast = noncoding_df.loc[noncoding_df['cancer'] != 'Breast']
        noncoding_nonbreast.to_csv(f"{path_file}noncoding_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        # Get the data frames for breast cancer and not breast cancer
        noncoding_breast = pd.read_csv(f"{path_file}noncoding_breast.tsv", sep='\t')
        noncoding_nonbreast = pd.read_csv(f"{path_file}noncoding_nonbreast.tsv", sep='\t')
    return noncoding_breast, noncoding_nonbreast


def filter_coding(filter_par, df_whole, path_file):
    """
    Filter or grab the coding data
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param df_whole: Dataframe with all the information out of the database
    :param path_file:  Path where the files are saved
    :return: coding_df: Dataframe with the coding information (SNPs)
    """
    if filter_par:
        # Get the coding SNPs
        coding_df = df_whole.loc[
            (df_whole['in_transcript'] == 1) | (df_whole['in_coding'] == 1) | (df_whole['in_exon'] == 1)]
        coding_df.to_csv(f"{path_file}coding_df.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        # Get the dataframe of the coding SNPs
        coding_df = pd.read_csv(f"{path_file}coding_df.tsv", sep='\t')
    return coding_df


def filter_coding_cancer(filter_par, coding_df, path_file):
    """
    Filter the coding data on breast cancer and nonbreast cancer
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param coding_df: Dataframe with the coding information (SNPs)
    :param path_file:  Path where the files are saved
    :return: coding_breast: All the coding SNPs of the breast cancer donors
             coding_nonbreast: All the coding SNPs of the nonbreast cancer donors
    """
    if filter_par:
        # Filters out breast cancer data and not breast cancer data
        coding_breast = coding_df.loc[coding_df['cancer'] == 'Breast']
        coding_breast.to_csv(f"{path_file}coding_breast.tsv", sep='\t', encoding='utf-8', index=False)
        coding_nonbreast = coding_df.loc[coding_df['cancer'] != 'Breast']
        coding_nonbreast.to_csv(f"{path_file}coding_nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    else:
        # Get the data frames for breast cancer and not breast cancer
        coding_breast = pd.read_csv(f"{path_file}coding_breast.tsv", sep='\t')
        coding_nonbreast = pd.read_csv(f"{path_file}coding_nonbreast.tsv", sep='\t')
    return coding_breast, coding_nonbreast


def get_all_data(filter_par, path_file, path_db):
    """
    Get the breast cancer and not breast cancer data
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param path_file:  Path where the files are saved
    :param path_db: Path to the database
    :return: all_breast: All the SNPs of the breast cancer donors
             all_nonbreast: All the SNPs of the nonbreast cancer donors
             num_donor_b: number of donors for breast cancer
             num_donor_nb: number of donors for nonbreast cancer
             len(all_breast): length of dataframe (breast)
             len(all_nonbreast): length of dataframe (nonbrast) 
    """
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        all_breast, all_nonbreast = filter_whole_cancer(filter_par, df_whole, path_file)
    else:
        all_breast, all_nonbreast = filter_whole_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(all_breast['donor_ID'])))
    num_donor_nb = len(list(set(all_nonbreast['donor_ID'])))

    return all_breast, all_nonbreast, num_donor_b, num_donor_nb, len(all_breast), len(all_nonbreast)


def get_noncoding_data(filter_par, path_file, path_db):
    """
    Get the breast cancer and not breast cancer data (noncoding SNPs)
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param path_file:  Path where the files are saved
    :param path_db: Path to the database
    :return: noncoding_breast: All the noncoding SNPs of the breast cancer donors
             noncoding_nonbreast All the noncoding SNPs of the nonbreast cancer donors
             num_donor_b: number of donors for breast cancer
             num_donor_nb: number of donors for nonbreast cancer
             len(noncoding_breast): length of dataframe (breast)
             len(noncoding_nonbreast): length of dataframe (nonbrast) 
    """
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        noncoding_df = filter_noncoding(filter_par, df_whole, path_file)
        noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(filter_par, noncoding_df, path_file)
    else:
        noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(noncoding_breast['donor_ID'])))
    num_donor_nb = len(list(set(noncoding_nonbreast['donor_ID'])))

    return noncoding_breast, noncoding_nonbreast, num_donor_b, num_donor_nb, len(noncoding_breast), len(
        noncoding_nonbreast)


def get_coding_data(filter_par, path_file, path_db):
    """
    Get the breast cancer and not breast cancer data (coding SNPs)
    :param filter_par: This parameter is true if it needs to retrieve data from the database and false if 
                       it has already done so and this data is stored somewhere. If false, it therefore uses
                       the previously stored data.
    :param path_file:  Path where the files are saved
    :param path_db: Path to the database
    :return: coding_breast: All the coding SNPs of the breast cancer donors
             coding_nonbreast: All the coding SNPs of the nonbreast cancer donors
             num_donor_b: number of donors for breast cancer
             num_donor_nb: number of donors for nonbreast cancer
             len(coding_breast): length of dataframe (breast)
             len(coding_nonbreast): length of dataframe (nonbrast) 
    """
    if filter_par:
        df_whole = get_data_db(filter_par, path_file, path_db)
        coding_df = filter_coding(filter_par, df_whole, path_file)
        coding_breast, coding_nonbreast = filter_coding_cancer(filter_par, coding_df, path_file)
    else:
        coding_breast, coding_nonbreast = filter_coding_cancer(filter_par, '', path_file)
    num_donor_b = len(list(set(coding_breast['donor_ID'])))
    num_donor_nb = len(list(set(coding_nonbreast['donor_ID'])))

    return coding_breast, coding_nonbreast, num_donor_b, num_donor_nb, len(coding_breast), len(coding_nonbreast)


def main():
    # Call get_config
    config = get_config('gearshift')
    path_file = config['analyse']
    path_db = config['database_get_data']
    # Get dataframe out of database
    df_whole = get_data_db(True, path_file, path_db)
    all_breast, all_nonbreast = filter_whole_cancer(True, df_whole, path_file)
    noncoding_df = filter_noncoding(True, df_whole, path_file)
    noncoding_breast, noncoding_nonbreast = filter_noncoding_cancer(True, noncoding_df, path_file)
    coding_df = filter_coding(True, df_whole, path_file)
    coding_breast, coding_nonbreast = filter_coding_cancer(True, coding_df, path_file)


if __name__ == '__main__':
    main()
