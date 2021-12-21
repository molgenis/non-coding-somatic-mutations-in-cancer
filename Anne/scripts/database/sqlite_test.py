import sqlite3
import glob
import pandas as pd


def create_db(cursor):
    """

    :param cursor:
    :return:
    """
    # Create tables
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `project`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `project`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `project_ID` VARCHAR(45) NOT NULL  
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `donor`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `donor_ID` VARCHAR(45) NOT NULL,
        `project_ID` INT NOT NULL,
        CONSTRAINT `fk_donor_project`
            FOREIGN KEY (`project_ID`)
            REFERENCES `project` (`ID`)
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `snp`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `chr` VARCHAR(45) NULL DEFAULT NULL,
        `pos_start` INT NULL DEFAULT NULL,
        `pos_end` INT NULL DEFAULT NULL,
        `ref` VARCHAR(1000) NULL DEFAULT NULL,
        `alt` VARCHAR(1000) NULL DEFAULT NULL,
        `genome_version` VARCHAR(45) NULL DEFAULT NULL,
        `depth` INT NULL DEFAULT NULL,
        `GLOC` TEXT CHECK( `GLOC` IN ('int','ext','non') ) NULL DEFAULT NULL,
        `platform` VARCHAR(45) NULL DEFAULT NULL,
        `seq_strategy` VARCHAR(45) NULL DEFAULT NULL,
        `tissue_id` VARCHAR(45) NULL DEFAULT NULL,
        UNIQUE (`chr`, `pos_start`, `pos_end`, `ref`, `alt`, `genome_version`, `GLOC`, `platform`, `seq_strategy`) 
        ON CONFLICT REPLACE
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `donor_has_snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor_has_snp`(
        `donor_ID` INT NOT NULL,
        `donor_project_ID` INT NOT NULL,
        `snp_ID` INT NOT NULL,
        CONSTRAINT `fk_donor_has_snp_donor1`
            FOREIGN KEY (`donor_ID`, `donor_project_ID`)
            REFERENCES `donor` (`ID`, `project_ID`)
        CONSTRAINT `fk_donor_has_snp_snp1`
            FOREIGN KEY (`snp_ID`)
            REFERENCES `snp` (`ID`)
    )
    """)
    return cursor


def fill_database(mydb_connection, cursor, df):
    """

    :param mydb_connection:
    :param cursor:
    :param df:
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        cursor.execute("""INSERT INTO project (project_ID) 
                          VALUES ('%s')""" % (str(project_id)))
        # Get the last ID (private key of the project table) used
        last_id_project = cursor.lastrowid
        print(f"{last_id_project} - {project_id}")
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        print(f"donors: {len(set(select_project['donor_id']))}")
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in list(set(select_project['donor_id'])):
            cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                              VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
            # Get the last ID (private key of the donor table) used
            last_id_donor = cursor.lastrowid
            # Filter dataframe on donor_id
            select_donor = select_project[select_project['donor_id'] == donor_id]
            # Loop over rows in dataframe (select_donor)
            for index, row in select_donor.iterrows():
                # See if an SNP already exists with these values
                cursor.execute(
                    """SELECT *
                    FROM snp
                    WHERE chr = '%s' AND pos_start = %s AND pos_end = %s AND 
                    ref = '%s' AND alt = '%s' AND genome_version = '%s' AND GLOC = '%s' 
                    AND platform = '%s' AND seq_strategy = '%s';""" %
                    (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                     str(row['ref']), str(row['alt']), str(row['genome_version']), str('int'),
                     str(row['platform']), str(row['seq_strategy'])))
                check_snp = cursor.fetchall()
                # If the SNP does not exist add it to the database
                if not check_snp:
                    cursor.execute("""
                        INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                                     GLOC, platform, seq_strategy, tissue_id)
                        VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s', '%s')""" %
                                   (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                                    str(row['ref']), str(row['alt']), str(row['genome_version']), int(row['depth']),
                                    str('int'), str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
                    # Get the last ID (private ket of the snp table) used
                    last_id_snp = cursor.lastrowid
                    # Fill the table donor_has_snp
                    cursor.execute("""
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
                        cursor.execute(
                            """SELECT *
                            FROM donor_has_snp
                            WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s;""" %
                            (int(last_id_project), int(last_id_donor), int(id_snp))
                        )
                        check_donor_snp = cursor.fetchall()
                        # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                        # with this combination
                        if not check_donor_snp:
                            cursor.execute("""
                                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                                VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(id_snp)))
    # Committing the current transactions
    mydb_connection.commit()


def read_files(path, mydb_connection, cursor):
    """

    :param path:
    :param mydb_connection:
    :param cursor:
    :return:
    """
    path_files = f'{path}*.tsv'
    # Loop over all files
    for fname in glob.glob(path_files):
        print(fname)
        # Read file
        df = pd.read_csv(fname, sep='\t')
        # Drop all duplicates (only depth and tissue_id may differ from all columns)
        df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
        fill_database(mydb_connection, cursor, df)


def main():
    """

    :return:
    """
    try:
        # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
        mydb_connection = sqlite3.connect(f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship.db')
        # Setting row_factory property of connection object to
        # sqlite3.Row(sqlite3.Row is an implementation of row_factory)
        mydb_connection.row_factory = sqlite3.Row
        # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
        cursor = mydb_connection.cursor()
        cursor = create_db(cursor)
        read_files('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/', mydb_connection, cursor)

        # cursor.execute('SELECT * FROM donor_has_snp')
        # for x in cursor:
        #     print(f"{x['donor_project_ID']} - {x['donor_ID']} - {x['snp_ID']}")
    except sqlite3.Error as er:
        print("Error while connecting to sqlite", er)
    finally:
        if mydb_connection:
            mydb_connection.close()
            print("The SQLite connection is closed")


main()
