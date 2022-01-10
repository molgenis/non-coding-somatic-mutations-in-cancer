import sqlite3
import pandas as pd
import sys
import io
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp

from db_ob import Database



# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    """

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def filter_add(db_path, df, alter):
    db = Database(db_path) #sys.argv[1]
    mydb_connection = db.mydb_connection
    cursor = db.cursor
    dbSNP = df[df['ID'].str.contains("rs")]
    print(dbSNP.head())
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
    

    path = "D:/Hanze_Groningen/STAGE/db/bdsnp filter/chr1_ann.vcf"
    # Read vcf file
    df = read_vcf(sys.argv[2])#(sys.argv[1].strip())
    df_shuffled = df.sample(frac=1)
    df_splits = np.array_split(df_shuffled, mp.cpu_count())
    arg_multi_list = []
    for df_s in df_splits:
        arg_multi_list.append((sys.argv[1], df_s, sys.argv[3]))

    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=filter_add, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # filter_add(df, mydb_connection, cursor, sys.argv[3])
    # db.count_values('ID_dbSNP', 'snp')
    # cursor.execute('SELECT * FROM snp')
    # for x in cursor:
    #     print(f"{x['ID_dbSNP']}")

    # db.close()

main()
