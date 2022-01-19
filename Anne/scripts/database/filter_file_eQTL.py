#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database


def fill_eQTL(eQTL_df, cursor, mydb_connection, alter):

    if alter == 'ALTER':
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `ID_eQTL` VARCHAR(45) NULL DEFAULT NULL
                        """)
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `eQTL` BOOLEAN DEFAULT(FALSE)
                        """)

    for index, row in eQTL_df.iterrows():
        cursor.execute(
                """UPDATE snp
                    SET ID_eQTL = '%s', eQTL = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    AND ref = '%s' AND alt = '%s';""" %
                (str(row['SNP']), str(row['SNPChr']), int(row['SNPPos']), 
                int(row['SNPPos']), str(row['AssessedAllele']), str(row['OtherAllele'])))
    mydb_connection.commit()


def main():
    db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    eQTL_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    db = Database(db_path) #sys.argv[1]
    eQTL_df = pd.read_csv(eQTL_path, sep='\t')
    fill_eQTL(eQTL_df, db.cursor, db.mydb_connection, 'ALTER')
    db.close()
    



if __name__ == '__main__':
    main()

