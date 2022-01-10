import sqlite3
import glob
import pandas as pd
import sys
import io
import os
import numpy as np
from multiprocessing import Pool
import multiprocessing as mp

from db_ob import Database


def check_gene(gene_df, path_db):
    """

    :param path_fgene:
    :return:
    """
    print('check_gene')
    db = Database(sys.argv[1]) #sys.argv[1]
    mydb_connection = db.mydb_connection
    cursor = db.cursor
    
    for index, row in gene_df.iterrows():
        print(index)
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
        print(f"COUNT - {row['exonCount']}")
        for i in range(int(row['exonCount'])):
            cursor.execute(
                """UPDATE snp
                    SET in_exon = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
            mydb_connection.commit()

    db.close()





def main():
    # path_fgene = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/snp132_ucsc_hg19_checkGene.bed'
    # path_db = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long2.db'
    # db = Database(sys.argv[1]) #sys.argv[1]
    # mydb_connection = db.mydb_connection
    # cursor = db.cursor

    gene_df = pd.read_csv(sys.argv[2], sep='\t')
    gene_shuffled = gene_df.sample(frac=1)
    gene_splits = np.array_split(gene_shuffled, mp.cpu_count())
    arg_gene = []
    for df_s in gene_splits:
        arg_gene.append((df_s, sys.argv[1]))
    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=check_gene, iterable=arg_gene)
    pool.close()
    pool.join()

    # db.close()

main()