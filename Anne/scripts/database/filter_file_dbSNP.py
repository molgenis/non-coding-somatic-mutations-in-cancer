import sqlite3
import pandas as pd
import sys
import io

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

def filter_add(df, mydb_connection, cursor):
    dbSNP = df[df['ID'].str.contains("rs")]
    print(dbSNP.head())
    # cursor.execute(f"""
    #                 ALTER TABLE snp
    #                 ADD `ID_dbSNP` VARCHAR(45) NULL DEFAULT NULL
    #                 """)
    print(len(dbSNP))
    for index, row in dbSNP.iterrows():
        # print(index)
        cursor.execute(
                """UPDATE snp
                    SET ID_dbSNP = '%s'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    AND ref = '%s' AND alt = '%s';""" %
                (str(row['ID']), str(row['CHROM'].replace('chr', '')), int(row['POS']), 
                int(row['POS']), str(row['REF']), str(row['ALT'])))
        mydb_connection.commit()


#     ALTER TABLE table_name  
#  ADD new_column_name column_definition 




def main():
    db = Database() #sys.argv[1]
    mydb_connection = db.mydb_connection
    cursor = db.cursor

    path = "D:/Hanze_Groningen/STAGE/db/bdsnp filter/chr1_ann.vcf"
    # Read vcf file
    df = read_vcf(path)#(sys.argv[1].strip())
    filter_add(df, mydb_connection, cursor)
    db.count_values('snp', 'ID_dbSNP')
    # cursor.execute('SELECT * FROM snp')
    # for x in cursor:
    #     print(f"{x['ID_dbSNP']}")

    db.close()

main()