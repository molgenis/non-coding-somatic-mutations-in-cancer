#!/usr/bin/env python3
from os import sep
import pandas as pd
import sys
import math # nan
import numpy as np

from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



def fill_project(db, df, project_cancer):
    # Project
    print('projects')
    project_id_set = list(set(df['project_id']))
    dict_project = dict()
    for project_id in project_id_set:
        cancer = project_cancer[project_id]
        db.cursor.execute("""INSERT INTO project (project_ID, cancer) 
                            VALUES ('%s', '%s')""" % (str(project_id), str(cancer)))
        last_id_project = db.cursor.lastrowid
        dict_project[project_id] = last_id_project
    return dict_project

def fill_donor(db, df, dict_project, donor_info):
    print('donors')
    donor_id_set = list(set(df['donor_id']))
    dict_donor = dict()
    for donor_id in donor_id_set:
        print(donor_id)
        select_items = df[df['donor_id'] == donor_id]
        project_id = list(set(select_items['project_id']))[0]
        db.cursor.execute("""INSERT INTO donor (donor_ID, project_ID, sex, vital_status, age_at_diagnosis, 
                                                    age_at_last_followup, disease_status_last_followup)
                                VALUES ('%s', %s, '%s', '%s', %s, %s, '%s')""" %
                            (str(donor_id), int(dict_project[project_id]),donor_info['donor_sex'][donor_id],
                            donor_info['donor_vital_status'][donor_id],
                            donor_info['donor_age_at_diagnosis'][donor_id],
                            donor_info['donor_age_at_last_followup'][donor_id],
                            donor_info['disease_status_last_followup'][donor_id]))
                
        # Get the last ID (private key of the donor table) used
        last_id_donor = db.cursor.lastrowid
        dict_donor[donor_id] = last_id_donor
    return dict_donor

def fill_snps(db, df):
    print('snps')
    df['snp'] = df['chr'] + '_' + df['pos_start'].map(str) + '_' + df['pos_end'].map(str) + '_' + df['ref'] + '_' + df['alt'] + '_' + df['genome_version'] + '_' + df['platform'] + '_' + df['seq_strategy']
    snp_set = list(set(df['snp']))
    dict_snp = dict()
    for snp in snp_set:
        snp_split = snp.split('_')
        db.cursor.execute("""
            INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version,
                        platform, seq_strategy)
            VALUES ('%s', %s, %s, '%s', '%s', '%s', '%s', '%s')""" %
                            (str(snp_split[0]), int(snp_split[1]), int(snp_split[2]),
                            str(snp_split[3]), str(snp_split[4]), str(snp_split[5]),
                            str(snp_split[6]), str(snp_split[7])))
        last_id_snp = db.cursor.lastrowid
        dict_snp[snp] = last_id_snp
    return dict_snp, df

def fill_donor_has_snp(db, df, dict_snp, dict_project, dict_donor, specimen_df):
    print('connection')
    df['ALL'] = df['project_id'] + '___' + df['donor_id'] + '___' + df['snp'] + '___' + df['specimen_id'] + '___' + df['total_read_count'].map(str) + '___' + df['mutant_allele_read_count'].map(str)
    all_set = list(set(df['ALL']))
    for all in all_set:
        all_split = all.split('___')
        # TODO # TODO # TODO # TODO
        specimen_type = specimen_df[specimen_df['icgc_specimen_id'] == all_split[3]]['specimen_type'].values[0]
        db.cursor.execute(
            """SELECT *
            FROM tissue
            WHERE specimen_type = '%s';""" %
            (str(specimen_type)))
        check_specimen = db.cursor.fetchall()
        for spe in check_specimen:
            # Get ID of the tissue
            tissue_id = int(spe['ID'])
        # TODO # TODO # TODO # TODO
        if all_split[4] == 'NULL':
            total_read = 'NULL'
        else:
            total_read = int(float(all_split[4]))

        if all_split[5] == 'NULL':
            mut_read = 'NULL'
        else:
            mut_read = int(float(all_split[5]))
            
        db.cursor.execute("""
                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID, tissue_id, specimen_id, 
                                            total_read_count, mutant_allele_read_count)
                VALUES (%s, %s, %s, %s, '%s', %s, %s)""" %
                                (int(dict_project[all_split[0]]), int(dict_donor[all_split[1]]), int(dict_snp[all_split[2]]), int(tissue_id),
                                str(all_split[3]), total_read, mut_read))


def read_file(path, db, project_cancer, donor_info, specimen_df):
    """
    Read the file and calls fill_database
    :param path: Path of the file to be read
    :param db: The database object
    :param project_cancer: Dictionary with as key project_id and as value kind of cancer
    :param donor_info: Dictionary with info about the donors. Key is column name and
                       value a dictionary with donor_ids and there info of that column.
    :param specimen_df: File with specimen_type (Normal or tumor tissue)
    :return: 
    """
    # Read file
    print('READ')
    df = pd.read_csv(path, sep='\t', compression='gzip')
    # Drop all duplicates
    print('remove duplicates')
    df = df.drop_duplicates()
    # Selects only the SNPs found with WGS
    print('filter WGS')
    df = df.loc[df['seq_strategy'] == 'WGS']
    print('fill')
    df[['total_read_count', 'mutant_allele_read_count']] = df[['total_read_count', 'mutant_allele_read_count']].fillna('NULL')
    # Checked if df is not empty (Some projects do not contain WGS SNPs)
    print('check')
    if len(df) != 0:
        # SNPs
        print('snps')
        dict_snp, df = fill_snps(db, df)
        db.mydb_connection.commit()
        # Project
        print('project')
        dict_project = fill_project(db, df, project_cancer)  
        db.mydb_connection.commit()      
        # Donor
        print('donor')
        dict_donor = fill_donor(db, df, dict_project, donor_info)
        db.mydb_connection.commit()
        # ALL
        print('all')
        fill_donor_has_snp(db, df, dict_snp, dict_project, dict_donor, specimen_df)        
    db.mydb_connection.commit()


def fill_tissue(specimen_df, db):
    """
    The tissue table is being filled.
    :param specimen_df: File with specimen_type (Normal or tumor tissue)
    :param db: The database object
    :return:
    """
    # Loop over set of specimen_type and add it to the database
    for specimen_type in list(set(specimen_df['specimen_type'])):
        db.cursor.execute(
            """SELECT *
                    FROM tissue
                    WHERE specimen_type = '%s';""" %
            (str(specimen_type)))
        check_tissue = db.cursor.fetchall()
        # Checked if tissue already exist in the database.
        # If it doesn't exist yet it will be added.
        if not check_tissue:
            # Check if tissue is normal tissue or tumor tissue.
            if 'normal' in specimen_type.lower():
                hc_tumor = 'FALSE'
            else:
                hc_tumor = 'TRUE'
            db.cursor.execute("""INSERT INTO tissue (specimen_type, type) 
                            VALUES ('%s', '%s')""" % (str(specimen_type), str(hc_tumor)))
    # Committing the current transactions
    db.mydb_connection.commit()

def main():
    config = get_config()
    # Path to file with project ID and kind of cancer
    site_path = config['site']  # "E:/STAGE/Site/site.csv"
    site_df = pd.read_csv(site_path, sep=';')
    # Make dictionary of site_df, and get only the column cancer
    project_cancer = site_df.set_index('project_ID').to_dict('dict')['cancer']
    # Path to file with donor information like age, sex etc.
    donor_info_path = config['donor_info']  # "E:/STAGE/WGS/donor.tsv"
    donor_info_df = pd.read_csv(donor_info_path, sep='\t')
    # Make of gender (sex) and vital status a boolean
    donor_info_df['donor_sex'] = donor_info_df['donor_sex'].map({'male': 'FALSE', 'female': 'TRUE'})
    donor_info_df['donor_vital_status'] = donor_info_df['donor_vital_status'].map(
        {'deceased': 'FALSE', 'alive': 'TRUE'})
    # Replace nan and empty values with NULL
    donor_info_df = donor_info_df.replace(np.nan,'NULL')
    donor_info_df = donor_info_df.replace('','NULL')
    # Make dictionary of donor_info_df
    donor_info = donor_info_df.set_index('icgc_donor_id').to_dict('dict')
    # Path to file with specimen_type (Normal or tumor tissue)
    specimen_path = config['all_specimen']  # "E:/STAGE/WGS/all_specimen.tsv"
    specimen_df = pd.read_csv(specimen_path, sep='\t')
    # Make Database object
    db = Database(config['test_database'])
    # Database connection
    # db = Database(path_db)
    # Calls fill_tissue
    print('fill_tissue')
    fill_tissue(specimen_df, db)
    # Calls read_file
    print('read_file')
    read_file(config['ALL_WGS_tot_mut'], db, project_cancer, donor_info, specimen_df)
    # read_file('D:\Hanze_Groningen\STAGE\DATAB\SKCM-US_db_NEW.tsv.gz', db) #, project_cancer, donor_info, specimen_df)

    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
