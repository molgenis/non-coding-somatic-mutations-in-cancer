#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp

from Database import Database


def check_gene(gene_df, mydb_connection, cursor):
    """

    :param path_fgene:
    :return:
    """
    print('check_gene')
    # Loop over gene_df
    for index, row in gene_df.iterrows():
        print(index)
        # Update in_transcript
        cursor.execute(
            """UPDATE snp
                SET in_transcript = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
        # mydb_connection.commit()
        # Update in_coding
        cursor.execute(
            """UPDATE snp
                SET in_coding = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
        # mydb_connection.commit()
        # Get start and end of the exons
        exon_start = row['exonStarts'].rstrip(',').split(',')
        exon_end = row['exonEnds'].rstrip(',').split(',')
        # print(f"COUNT - {row['exonCount']}")
        # Loop over the exons start-end
        for i in range(int(row['exonCount'])):
            # Update in_exon
            cursor.execute(
                """UPDATE snp
                    SET in_exon = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
            # mydb_connection.commit()
    # Add to database
    mydb_connection.commit()


def main():
    # path_fgene = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/snp132_ucsc_hg19_checkGene.bed'
    # path_db = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long2.db'
    # db = Database(sys.argv[1]) #sys.argv[1]
    # mydb_connection = db.mydb_connection
    # cursor = db.cursor
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db"
    # Path of the genes and there positions
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')  # sys.argv[2], sep='\t')
    # Database connection
    db = Database(path_db)  # sys.argv[1]
    # Call check_gene
    check_gene(gene_df, db.mydb_connection, db.cursor)
    # Close connections
    db.close()


if __name__ == '__main__':
    main()
