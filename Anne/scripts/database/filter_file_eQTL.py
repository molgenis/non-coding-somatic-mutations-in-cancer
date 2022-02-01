#!/usr/bin/env python3
import pandas as pd
import sys
from Database import Database


def fill_eQTL_M(eQTL_df, cursor, mydb_connection, alter, ID_eQT, eQT):
    """
    """
    # When alter is ALTER, it means that nothing has been changed in the database and that a column has to be added to the database.
    if alter == 'ALTER':
        cursor.execute("""
                        ALTER TABLE snp
                        ADD `%s` VARCHAR(45) NULL DEFAULT NULL
                        """ %
                        (ID_eQT))
        cursor.execute("""
                        ALTER TABLE snp
                        ADD `%s` BOOLEAN DEFAULT(FALSE)
                        """ %
                        (eQT))
    # Loop over rows in eQTL file
    for index, row in eQTL_df.iterrows():
        # Add ID_eQTL and eQTL to snp that matches the matches in chr, pos_start, pos_end, ref and alt.
        cursor.execute(
                """UPDATE snp
                    SET %s = '%s', %s = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    AND ref = '%s' AND alt = '%s';""" %
                (ID_eQT, eQT, str(row['SNP']), str(row['SNPChr']), int(row['SNPPos']), 
                int(row['SNPPos']), str(row['AssessedAllele']), str(row['OtherAllele'])))
    # Add to database
    mydb_connection.commit()


def main():
    # Database file
    db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    # File with eQTLs
    eQTL_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # eQTM_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # Database class
    db = Database(db_path) #sys.argv[1]
    # Read file
    eQTL_df = pd.read_csv(eQTL_path, sep='\t')
    # eQTM_df = pd.read_csv(eQTL_path, sep='\t')
    # Call fill_eQTL
    fill_eQTL_M(eQTL_df, db.cursor, db.mydb_connection, 'ALTER', 'ID_eQTL', 'eQTL')
    # fill_eQTL_M(eQTL_df, db.cursor, db.mydb_connection, 'ALTER', 'ID_eQTM', 'eQTM')
    # Close connection cursor and mydb_connection
    db.close()


if __name__ == '__main__':
    main()

