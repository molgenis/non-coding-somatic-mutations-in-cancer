#!/usr/bin/env python3
import pandas as pd
import sys

from Database import Database


def fill_database(df, db):
    """
    Fill the database with the info in the df
    :param df:  Dataframe of a project
    :param db:  the database object
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        print(project_id)
        # Fill project table
        db.cursor.execute("""INSERT INTO project (project_ID) 
                        VALUES ('%s')""" % (str(project_id)))
        # Get the last ID (private key of the project table) used
        last_id_project = db.cursor.lastrowid
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in list(set(select_project['donor_id'])):
            # Fill donor table
            db.cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                            VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
            # Get the last ID (private key of the donor table) used
            last_id_donor = db.cursor.lastrowid
            # Filter dataframe on donor_id
            select_donor = select_project[select_project['donor_id'] == donor_id]
            # Loop over rows in dataframe (select_donor)
            for index, row in select_donor.iterrows():
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
                        INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                                    platform, seq_strategy, tissue_id)
                        VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s')""" %
                                      (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                                       str(row['ref']), str(row['alt']), str(row['genome_version']),
                                       int(row['depth']),
                                       str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
                    # Get the last ID (private ket of the snp table) used
                    last_id_snp = db.cursor.lastrowid
                    # Fill donor_has_snp table
                    db.cursor.execute("""
                        INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                        VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(last_id_snp)))
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
                            WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s;""" %
                            (int(last_id_project), int(last_id_donor), int(id_snp))
                        )
                        check_donor_snp = db.cursor.fetchall()
                        # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                        # with this combination
                        if not check_donor_snp:
                            # Fill donor_has_snp table
                            db.cursor.execute("""
                                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                                VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(id_snp)))
    # Committing the current transactions
    db.mydb_connection.commit()


def read_file(path, db):
    """
    Read the file and calls fill_database
    :param path: Path of the file to be read
    :param db: the database object
    :return: 
    """
    # Read file
    df = pd.read_csv(path, sep='\t')
    # Drop all duplicates (only depth and tissue_id may differ from all columns)
    df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
    # Calls fill_database
    fill_database(df, db)


def main():
    # Make Database object
    db = Database(sys.argv[1])
    # Calls read_file
    read_file(sys.argv[2], db)
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
