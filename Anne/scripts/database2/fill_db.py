#!/usr/bin/env python3
from os import sep
import pandas as pd
import sys
import math # nan
import numpy as np

from Database import Database


def fill_snp_tissue_donorsnp(db, select_donor, specimen_df, last_id_project, last_id_donor):
    """
    Fill the database table `snp`, `donor_has_snp` with the info in the df
    :param db:  The database object
    :param select_donor:  Dataframe with only one (selected) donor
    :param specimen_df: File with specimen_type (Normal or tumor tissue)
    :param last_id_project: ID in the database of (selected) project
    :param last_id_donor: ID in the database of (selected) donor
    :return:
    """
    # Loop over rows in dataframe (select_donor)
    for index, row in select_donor.iterrows():
        # Get specimen_id
        specimen_id = row['specimen_id']
        # Get specimen_tye
        specimen_type = specimen_df[specimen_df['icgc_specimen_id'] == specimen_id]['specimen_type'].values[0]
        db.cursor.execute(
            """SELECT *
            FROM tissue
            WHERE specimen_type = '%s';""" %
            (str(specimen_type)))
        check_specimen = db.cursor.fetchall()
        for spe in check_specimen:
            # Get ID of the tissue
            tissue_id = int(spe['ID'])
        # See if an SNP already exists with these values
        db.cursor.execute(
            """SELECT *
            FROM snp
            WHERE chr = '%s' AND pos_start = %s AND pos_end = %s AND 
            ref = '%s' AND alt = '%s' AND genome_version = '%s' 
            AND platform = '%s' AND seq_strategy = '%s';""" %
            (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
             str(row['ref']), str(row['alt']), str(row['genome_version']),
             str(row['platform']), str(row['seq_strategy'])))
        check_snp = db.cursor.fetchall()
        # If the SNP does not exist add it to the database
        if not check_snp:
            # Fill snp table
            db.cursor.execute("""
                INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version,
                            platform, seq_strategy)
                VALUES ('%s', %s, %s, '%s', '%s', '%s', '%s', '%s')""" %
                              (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                               str(row['ref']), str(row['alt']), str(row['genome_version']),
                               str(row['platform']), str(row['seq_strategy'])))
            # Get the last ID (private ket of the snp table) used
            last_id_snp = db.cursor.lastrowid
            # Fill donor_has_snp table
            db.cursor.execute("""
                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID, tissue_id, specimen_id, total_read_count, mutant_allele_read_count)
                VALUES (%s, %s, %s, %s, '%s', %s, %s)""" %
                              (int(last_id_project), int(last_id_donor), int(last_id_snp), int(tissue_id),
                               str(specimen_id), row['total_read_count'], row['mutant_allele_read_count']))
        # If the snp already exists insert the link between the donor and the snp by filling in
        # the donor_has_snp table
        else:
            # Loop over the snp(s) it corresponds to, if all is well this is always 1 snp.
            for info in check_snp:
                # Get ID of the snp
                id_snp = int(info['ID'])
                # Check whether the combination donor and snp already exists
                db.cursor.execute(
                    """SELECT *
                    FROM donor_has_snp
                    WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s AND tissue_id = %s AND 
                    specimen_id = '%s' AND total_read_count = %s AND mutant_allele_read_count = %s;""" %
                    (int(last_id_project), int(last_id_donor), int(id_snp), int(tissue_id), str(specimen_id), 
                    row['total_read_count'], row['mutant_allele_read_count'])
                )
                check_donor_snp = db.cursor.fetchall()
                # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                # with this combination
                if not check_donor_snp:
                    # Fill donor_has_snp table
                    db.cursor.execute("""
                        INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID, tissue_id, specimen_id, 
                                                    total_read_count, mutant_allele_read_count)
                        VALUES (%s, %s, %s, %s, '%s', %s, %s)""" %
                                      (int(last_id_project), int(last_id_donor), int(id_snp), int(tissue_id),
                                       str(specimen_id), row['total_read_count'], row['mutant_allele_read_count']))


def fill_donor(db, select_project, donor_info, last_id_project, specimen_df):
    """
    Fill the database table `donor` with the info in the df
    :param db: The database object
    :param select_project: Dataframe with only one (selected) project
    :param donor_info: Dictionary with info about the donors. Key is column name and
                       value a dictionary with donor_ids and there info of that column.
    :param last_id_project: ID in the database of (selected) project
    :param specimen_df:  File with specimen_type (Normal or tumor tissue)
    :return:
    """
    # Loop over set of donor_ids in (last) project_id and add it to the database
    for donor_id in list(set(select_project['donor_id'])):
        db.cursor.execute(
            """SELECT *
                    FROM donor
                    WHERE donor_ID = '%s';""" %
            (str(donor_id)))
        check_donor = db.cursor.fetchall()
        # If the donor does not exist add it to the database
        if not check_donor:
            # print(f"{donor_id} donor_sex - {donor_info['donor_sex'][donor_id]}")
            # print(f"{donor_id} donor_vital_status - {donor_info['donor_vital_status'][donor_id]}")
            # print(f"{donor_id} donor_age_at_diagnosis - {donor_info['donor_age_at_diagnosis'][donor_id]}")
            # print(f"{donor_id} donor_age_at_last_followup- {donor_info['donor_age_at_last_followup'][donor_id]}")
            # print(f"{donor_id} disease_status_last_followup- {donor_info['disease_status_last_followup'][donor_id]}")
            # print('-----------')
            # Fill donor table
            db.cursor.execute("""INSERT INTO donor (donor_ID, project_ID, sex, vital_status, age_at_diagnosis, 
                                                    age_at_last_followup, disease_status_last_followup)
                                VALUES ('%s', %s, '%s', '%s', %s, %s, '%s')""" %
                            (str(donor_id), int(last_id_project), donor_info['donor_sex'][donor_id],
                            donor_info['donor_vital_status'][donor_id],
                            donor_info['donor_age_at_diagnosis'][donor_id],
                            donor_info['donor_age_at_last_followup'][donor_id],
                            donor_info['disease_status_last_followup'][donor_id]))
                
            # Get the last ID (private key of the donor table) used
            last_id_donor = db.cursor.lastrowid
        else:
            for inf_don in check_donor:
                # Get ID of the donor
                last_id_donor = int(inf_don['ID'])
        # Filter dataframe on donor_id
        select_donor = select_project[select_project['donor_id'] == donor_id]
        # Calls fill_snp_tissue_donorsnp
        fill_snp_tissue_donorsnp(db, select_donor, specimen_df, last_id_project, last_id_donor)


def fill_database(df, db, project_cancer, donor_info, specimen_df):
    """
    Fill the database with the info in the df
    :param df:  Dataframe of a project
    :param db:  The database object
    :param project_cancer: Dictionary with as key project_id and as value kind of cancer
    :param donor_info: Dictionary with info about the donors. Key is column name and
                       value a dictionary with donor_ids and there info of that column.
    :param specimen_df: File with specimen_type (Normal or tumor tissue)
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        print(project_id)
        # Selected the kind of cancer of that project(_id)
        cancer = project_cancer[project_id]
        db.cursor.execute(
            """SELECT *
                    FROM project
                    WHERE project_ID = '%s';""" %
            (str(project_id)))
        check_project = db.cursor.fetchall()
        # If the project does not exist add it to the database
        if not check_project:
            # Fill project table
            db.cursor.execute("""INSERT INTO project (project_ID, cancer) 
                            VALUES ('%s', '%s')""" % (str(project_id), str(cancer)))
            # Get the last ID (private key of the project table) used
            last_id_project = db.cursor.lastrowid
        else:
            for inf_pro in check_project:
                # Get ID of the donor
                last_id_project = int(inf_pro['ID'])
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        # Calls fill_donor
        fill_donor(db, select_project, donor_info, last_id_project, specimen_df)
    # Committing the current transactions
    db.mydb_connection.commit()


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
    df = pd.read_csv(path, sep='\t', compression='gzip')
    # Drop all duplicates
    df = df.drop_duplicates()
    # Selects only the SNPs found with WGS
    df = df.loc[df['seq_strategy'] == 'WGS']
    # Checked if df is not empty (Some projects do not contain WGS SNPs)
    if len(df) != 0:
        # Calls fill_database
        fill_database(df, db, project_cancer, donor_info, specimen_df)


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
    # Path to file with project ID and kind of cancer
    site_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/site.csv'  # "E:/STAGE/Site/site.csv"
    site_df = pd.read_csv(site_path, sep=';')
    # Make dictionary of site_df, and get only the column cancer
    project_cancer = site_df.set_index('project_ID').to_dict('dict')['cancer']
    # Path to file with donor information like age, sex etc.
    donor_info_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/donor.tsv'  # "E:/STAGE/WGS/donor.tsv"
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
    specimen_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/all_specimen.tsv'  # "E:/STAGE/WGS/all_specimen.tsv"
    specimen_df = pd.read_csv(specimen_path, sep='\t')
    # Make Database object
    db = Database(sys.argv[1])
    # Calls fill_tissue
    fill_tissue(specimen_df, db)
    # Calls read_file
    read_file(sys.argv[2], db, project_cancer, donor_info, specimen_df)
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
