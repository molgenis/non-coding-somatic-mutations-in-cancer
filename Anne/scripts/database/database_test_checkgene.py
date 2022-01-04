import sqlite3
import glob
import pandas as pd
import sys


def check_gene(path_fgene, mydb_connection, cursor):
    """

    :param path_fgene:
    :param mydb_connection:
    :param cursor:
    :return:
    """
    print('check_gene')
    gene_df = pd.read_csv(path_fgene, sep='\t')
    for index, row in gene_df.iterrows():
        cursor.execute(
            """UPDATE snp
                SET in_transcript = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
        mydb_connection.commit()
        cursor.execute(
            """UPDATE snp
                SET in_coding = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
        mydb_connection.commit()
        exon_start = row['exonStarts'].rstrip(',').split(',')
        exon_end = row['exonEnds'].rstrip(',').split(',')
        for i in range(int(row['exonCount'])):
            cursor.execute(
                """UPDATE snp
                    SET in_exon = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
            mydb_connection.commit()


def main():
    """

    :return:
    """
    path_fgene = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/snp132_ucsc_hg19_checkGene.bed'
    path_db = f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/{sys.argv[1]}.db'
    try:
        # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
        mydb_connection = sqlite3.connect(path_db)
        # Setting row_factory property of connection object to
        # sqlite3.Row(sqlite3.Row is an implementation of row_factory)
        mydb_connection.row_factory = sqlite3.Row
        # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
        cursor = mydb_connection.cursor()
        
        check_gene(path_fgene, mydb_connection, cursor)

        cursor.execute('SELECT * FROM snp')
        for x in cursor:
            print(f"{x['in_transcript']} - {x['in_coding']} - {x['in_exon']}")
    except sqlite3.Error as er:
        print("Error while connecting to sqlite", er)
    finally:
        if mydb_connection:
            mydb_connection.close()
            print("The SQLite connection is closed")


main()
print(f'END {sys.argv[1]}')
