#!/usr/bin/env python3
import pandas as pd
import sys
from Database import Database


def add_value(db, ID_eQT, eQT, close_eQT):
    """

    """
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` VARCHAR(45) NULL DEFAULT NULL
                    """ %
                    (ID_eQT))
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` BOOLEAN DEFAULT(FALSE)
                    """ %
                    (eQT))
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` BOOLEAN DEFAULT(FALSE)
                    """ %
                    (close_eQT))

def set_value(db, ID_eQT, eQT, row):
    """

    """
    # Add ID_eQTL and eQTL to snp that matches the matches in chr, pos_start, pos_end, ref and alt.
    db.cursor.execute(
            """UPDATE snp
                SET %s = '%s', %s = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                AND ref = '%s' AND alt = '%s';""" %
            (ID_eQT, str(row['SNP']), eQT, str(row['SNPChr']), int(row['SNPPos']), 
            int(row['SNPPos']), str(row['AssessedAllele']), str(row['OtherAllele'])))


def set_value_close_to(db, row, close_eQT, region):
    """

    """
    db.cursor.execute(
            """UPDATE snp
                SET %s = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                AND ref = '%s' AND alt = '%s';""" %
            (close_eQT, str(row['SNPChr']), int(row['SNPPos'])-region, 
            int(row['SNPPos'])+region, str(row['AssessedAllele']), str(row['OtherAllele'])))


def loop_eQTL(db, eQTL_df, ID_eQT, eQT, close_eQT, region):
    """

    """
    for index, row in eQTL_df.iterrows():
        print(index)
        set_value(db, ID_eQT, eQT, row)
        set_value_close_to(db, row, close_eQT, region)
    # Add to database
    db.mydb_connection.commit()


def main():
    # Database file
    db_path="D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"#"D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    # File with eQTLs
    eQTL_path = "D:/Hanze_Groningen/STAGE/eQTL/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # eQTM_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # Database class
    db = Database(db_path) #sys.argv[1]
    # Read file
    eQTL_df = pd.read_csv(eQTL_path, sep='\t')
    # eQTM_df = pd.read_csv(eQTL_path, sep='\t')
    region = 100
    # add_value(db, 'ID_eQTL', 'eQTL', 'close_eQTL')
    # Add to database
    db.mydb_connection.commit()
    # Call fill_eQTL
    loop_eQTL(db, eQTL_df, 'ID_eQTL', 'eQTL', 'close_eQTL', region)
    # Add to database
    db.mydb_connection.commit()
    # Close connection cursor and mydb_connection
    db.close()


if __name__ == '__main__':
    main()

