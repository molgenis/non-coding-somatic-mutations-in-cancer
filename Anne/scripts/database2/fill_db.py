#!/usr/bin/env python3
from os import sep
import pandas as pd
import sys

from Database import Database


def process_info_donor(donor_info, donor_id):
    if donor_info['donor_sex'][donor_id] == 'male':
        sex = str('FALSE')
    elif donor_info['donor_sex'][donor_id] == 'female':
        sex = str('TRUE')
    else:
        sex = str('FALSE')

    if donor_info['donor_vital_status'][donor_id] == 'deceased':
        vital_status = str('FALSE')
    elif donor_info['donor_vital_status'][donor_id] == 'alive':
        vital_status = str('TRUE')
    else:
        vital_status = str('FALSE')

    return sex, vital_status


def fill_database(df, db, project_cancer, donor_info, specimen_df):
    """
    Fill the database with the info in the df
    :param df:  Dataframe of a project
    :param db:  the database object
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        print(project_id)
        cancer = project_cancer[project_id]
        db.cursor.execute(
                    """SELECT *
                    FROM project
                    WHERE project_ID = '%s';""" %
                    (str(project_id)))
        check_project = db.cursor.fetchall()
        # If the SNP does not exist add it to the database
        if not check_project:
            # print('IF PROJECT')
            # Fill project table
            db.cursor.execute("""INSERT INTO project (project_ID, cancer) 
                            VALUES ('%s', '%s')""" % (str(project_id), str(cancer)))
            # Get the last ID (private key of the project table) used
            last_id_project = db.cursor.lastrowid
        else:
            # print('ELSE PROJECT')
            for inf_pro in check_project:
                # Get ID of the snp
                last_id_project = int(inf_pro['ID'])
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in list(set(select_project['donor_id'])):
            sex, vital_status = process_info_donor(donor_info, donor_id)
            db.cursor.execute(
                        """SELECT *
                        FROM donor
                        WHERE donor_ID = '%s';""" %
                        (str(donor_id)))
            check_donor = db.cursor.fetchall()
            # If the SNP does not exist add it to the database
            if not check_donor:
                # print('IF DONOR')
                # Fill donor table
                db.cursor.execute("""INSERT INTO donor (donor_ID, project_ID, sex, vital_status, age_at_diagnosis, age_at_last_followup, disease_status_last_followup)
                                VALUES ('%s', %s, '%s', '%s', %s, %s, '%s')""" % 
                                (str(donor_id), int(last_id_project), str(sex), str(vital_status),
                                int(donor_info['donor_age_at_diagnosis'][donor_id]), int(donor_info['donor_age_at_last_followup'][donor_id]),
                                str(donor_info['disease_status_last_followup'][donor_id])))
                # Get the last ID (private key of the donor table) used
                last_id_donor = db.cursor.lastrowid
            else:
                # print('ELSE DONOR')
                for inf_don in check_donor:
                    # Get ID of the snp
                    last_id_donor = int(inf_don['ID'])
            # Filter dataframe on donor_id
            select_donor = select_project[select_project['donor_id'] == donor_id]
            # Loop over rows in dataframe (select_donor)
            for index, row in select_donor.iterrows():
                specimen_id = row['specimen_id']
                specimen_type = specimen_df[specimen_df['icgc_specimen_id'] == specimen_id]['specimen_type'].values[0]
                # print(specimen_type)
                db.cursor.execute(
                    """SELECT *
                    FROM tissue
                    WHERE specimen_type = '%s';""" %
                    (str(specimen_type)))
                check_specimen = db.cursor.fetchall()
                for spe in check_specimen:
                        # Get ID of the snp
                        tissue_id = int(spe['ID'])
                # print(tissue_id)



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
                        INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID, tissue_id, specimen_id)
                        VALUES (%s, %s, %s, %s, '%s')""" %
                        (int(last_id_project), int(last_id_donor), int(last_id_snp), int(tissue_id), str(specimen_id)))
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
                            WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s AND tissue_id = %s AND specimen_id = '%s';""" %
                            (int(last_id_project), int(last_id_donor), int(id_snp), int(tissue_id), str(specimen_id))
                        )
                        check_donor_snp = db.cursor.fetchall()
                        # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                        # with this combination
                        if not check_donor_snp:
                            # Fill donor_has_snp table
                            db.cursor.execute("""
                                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID, tissue_id, specimen_id)
                                VALUES (%s, %s, %s, %s, '%s')""" % 
                                (int(last_id_project), int(last_id_donor), int(id_snp), int(tissue_id), str(specimen_id)))
    # Committing the current transactions
    db.mydb_connection.commit()



def read_file(path, db, project_cancer, donor_info, specimen_df):
    """
    Read the file and calls fill_database
    :param path: Path of the file to be read
    :param db: the database object
    :return: 
    """
    # Read file
    df = pd.read_csv(path, sep='\t')
    # Drop all duplicates (only depth and tissue_id may differ from all columns)
    df = df.drop_duplicates() #subset=df.columns.difference(['depth', 'tissue_id'])
    df = df.loc[df['seq_strategy'] == 'WGS']

    if len(df) != 0:
        # Calls fill_database
        fill_database(df, db, project_cancer, donor_info, specimen_df) #df, db, project_cancer, donor_info, specimen_df

def fill_tissue(specimen_df, db):
    for specimen_type in list(set(specimen_df['specimen_type'])):
        db.cursor.execute(
                    """SELECT *
                    FROM tissue
                    WHERE specimen_type = '%s';""" %
                    (str(specimen_type)))
        check_tissue = db.cursor.fetchall()
        # If the SNP does not exist add it to the database
        if not check_tissue:
            print('ADD')
            if 'normal' in specimen_type.lower():
                hc_tumor = 'FALSE'
            else:
                hc_tumor = 'TRUE'
            db.cursor.execute("""INSERT INTO tissue (specimen_type, type) 
                            VALUES ('%s', '%s')""" % (str(specimen_type), str(hc_tumor)))
    db.mydb_connection.commit()


def main():
    site_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/site.csv' #"E:/STAGE/Site/site.csv"
    site_df = pd.read_csv(site_path, sep=';')
    project_cancer = site_df.set_index('project_ID').to_dict('dict')['cancer']

    donor_info_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/donor.tsv' #"E:/STAGE/WGS/donor.tsv"
    donor_info_df = pd.read_csv(donor_info_path, sep='\t')
    donor_info = donor_info_df.set_index('icgc_donor_id').to_dict('dict')

    specimen_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/all_specimen.tsv' #"E:/STAGE/WGS/all_specimen.tsv"
    specimen_df = pd.read_csv(specimen_path, sep='\t')
    #SP85685
    # print(specimen_df[specimen_df['icgc_specimen_id'] == 'SP85685'])
    # Make Database object
    db = Database(sys.argv[1])

    fill_tissue(specimen_df, db)  
    # db.cursor.execute(
    #                 """SELECT *
    #                 FROM tissue
    #                 WHERE specimen_type = '%s';""" %
    #                 (str('Normal - blood derived')))
    # check_specimen = db.cursor.fetchall()
    # for x in check_specimen:
    #     # Get ID of the snp
    #     print(int(x['ID']))  
    #     y = int(x['ID'])
    # print(f'y - {y}')
    # Calls read_file
    read_file(sys.argv[2], db, project_cancer, donor_info, specimen_df)

    # db.cursor.execute(
    #                 """SELECT *
    #                 FROM donor;""")
    # y = db.cursor.fetchall()
    # for x in y:
    #     # Get ID of the snp
    #     print(int(x['ID']))

    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
