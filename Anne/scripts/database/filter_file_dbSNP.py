#!/usr/bin/env python3
import pandas as pd
import sys
import numpy as np
from Database import Database

# Also takes the folder 1 higher, so that I can do the import after
# sys.path.append("..")
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf


def filter_add(mydb_connection, cursor, df, alter):
    """
    """
    # When alter is ALTER, it means that nothing has been changed in the database and that a column has to be added to the database.
    if alter == 'ALTER':
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `ID_dbSNP` VARCHAR(45) NULL DEFAULT NULL
                        """)
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `germline` BOOLEAN DEFAULT(FALSE)
                        """)
    
    # Only grab the snps that contain 'rs' in the ID column. These snps are then known in the dbSNP.
    dbSNP = df[df['ID'].str.contains("rs")]
    print(len(dbSNP))
    # Loop over rows in dbSNP file
    for index, row in dbSNP.iterrows():
        print(index)
        # Add germline and ID_dbsnp to snp that matches the matches in chr, pos_start, pos_end, ref and alt.
        cursor.execute(
                """UPDATE snp
                    SET ID_dbSNP = '%s', germline = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    AND ref = '%s' AND alt = '%s';""" %
                (str(row['ID']), str(row['CHROM'].replace('chr', '')), int(row['POS']), 
                int(row['POS']), str(row['REF']), str(row['ALT'])))
        print('hoi')
    # Add to database
    mydb_connection.commit()


def main():
    # File with dbSNPs
    path_ann = "D:/Hanze_Groningen/STAGE/db/bdsnp filter/chr15_ann - kopie.vcf"
    # Database file
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (2).db"
    # Read vcf file
    df = read_vcf(path_ann)
    # Database class
    db = Database(path_db)
    # Call filter_add
    filter_add(db.mydb_connection, db.cursor, df, '')
    # Close connection cursor and mydb_connection
    db.close()


if __name__ == '__main__':
    main()
