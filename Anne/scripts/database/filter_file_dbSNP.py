#!/usr/bin/env python3
import pandas as pd
import sys
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp

from Database import Database

# Also takes the folder 1 higher, so that I can do the import after
sys.path.append("..")
# sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf


def filter_add(mydb_connection, cursor, df, alter):
   
    dbSNP = df[df['ID'].str.contains("rs")]
    
    if alter == 'ALTER':
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `ID_dbSNP` VARCHAR(45) NULL DEFAULT NULL
                        """)
        cursor.execute(f"""
                        ALTER TABLE snp
                        ADD `germline` BOOLEAN DEFAULT(FALSE)
                        """)
    print(len(dbSNP))
    for index, row in dbSNP.iterrows():
        cursor.execute(
                """UPDATE snp
                    SET ID_dbSNP = '%s', germline = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    AND ref = '%s' AND alt = '%s';""" %
                (str(row['ID']), str(row['CHROM'].replace('chr', '')), int(row['POS']), 
                int(row['POS']), str(row['REF']), str(row['ALT'])))
    mydb_connection.commit()


def main():
    path_ann = "D:/Hanze_Groningen/STAGE/db/bdsnp filter/chr15_ann - kopie.vcf"
    path_db = "D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    # Read vcf file
    df = read_vcf(path_ann) #sys.argv[2])#(sys.argv[1].strip())
    db = Database(path_db)
    filter_add(db.mydb_connection, db.cursor, df, '')
    db.close()

if __name__ == '__main__':
    main()
