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


def fill_database(df, db):
    """

    :param df:
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        print(project_id)
        db.cursor.execute("""INSERT INTO project (project_ID) 
                        VALUES ('%s')""" % (str(project_id)))
        # Committing the current transactions
        db.mydb_connection.commit()
        # Get the last ID (private key of the project table) used
        last_id_project = db.cursor.lastrowid
        print(f"{last_id_project} - {project_id}")
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        print(f"donors: {len(set(select_project['donor_id']))}")
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in list(set(select_project['donor_id'])):
            db.cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                            VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
            # Committing the current transactions
            db.mydb_connection.commit()
            # Get the last ID (private key of the donor table) used
            last_id_donor = db.cursor.lastrowid
            # Filter dataframe on donor_id
            select_donor = select_project[select_project['donor_id'] == donor_id]
            # Loop over rows in dataframe (select_donor)
            for index, row in select_donor.iterrows():
                # See if an SNP already exists with these values
                db.cursor.execute(
                    """SELECT *
                    FROM snp
                    WHERE chr = '%s' AND pos_start = %s AND pos_end = %s AND 
                    ref = '%s' AND alt = '%s' AND genome_version = '%s' 
                    AND platform = '%s' AND seq_strategy = '%s';""" %
                    (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                        str(row['ref']), str(row['alt']), str(row['genome_version']),
                        str(row['platform']), str(row['seq_strategy'])))
                check_snp = db.cursor.fetchall()
                # If the SNP does not exist add it to the database
                if not check_snp:
                    db.cursor.execute("""
                        INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                                    platform, seq_strategy, tissue_id)
                        VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s')""" %
                                        (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                                            str(row['ref']), str(row['alt']), str(row['genome_version']),
                                            int(row['depth']),
                                            str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
                    # Get the last ID (private ket of the snp table) used
                    last_id_snp = db.cursor.lastrowid
                    # Fill the table donor_has_snp
                    db.cursor.execute("""
                        INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                        VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(last_id_snp)))
                    # Committing the current transactions
                    db.mydb_connection.commit()
                # If the snp already exists insert the link between the donor and the snp by filling in
                # the donor_has_snp table
                else:
                    # Loop over the snp(s) it corresponds to, if all is well this is always 1 snp.
                    for info in check_snp:
                        # Get ID of the snp
                        id_snp = int(info['ID'])
                        # Check whether the combination donor and snp already exists
                        db.cursor.execute(
                            """SELECT *
                            FROM donor_has_snp
                            WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s;""" %
                            (int(last_id_project), int(last_id_donor), int(id_snp))
                        )
                        check_donor_snp = db.cursor.fetchall()
                        # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                        # with this combination
                        if not check_donor_snp:
                            db.cursor.execute("""
                                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                                VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(id_snp)))
                        # Committing the current transactions
                        db.mydb_connection.commit()

def check_gene(gene_df, db):
    """

    :param path_fgene:
    :return:
    """
    print('check_gene')
    
    for index, row in gene_df.iterrows():
        print(index)
        db.cursor.execute(
            """UPDATE snp
                SET in_transcript = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
        db.mydb_connection.commit()
        db.cursor.execute(
            """UPDATE snp
                SET in_coding = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
        db.mydb_connection.commit()
        exon_start = row['exonStarts'].rstrip(',').split(',')
        exon_end = row['exonEnds'].rstrip(',').split(',')
        print(f"COUNT - {row['exonCount']}")
        for i in range(int(row['exonCount'])):
            db.cursor.execute(
                """UPDATE snp
                    SET in_exon = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
            db.mydb_connection.commit()


def read_file(path, db):
    """
    
    :param path: 
    :param db: 
    :return: 
    """
    df = pd.read_csv(path, sep='\t')
    # Drop all duplicates (only depth and tissue_id may differ from all columns)
    df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
    df_shuffled = df.sample(frac=1)
    df_splits = np.array_split(df_shuffled, mp.cpu_count())
    arg_multi_list = []
    for df_s in df_splits:
        arg_multi_list.append((df_s, db))

    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=fill_database, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # fill_database(df, db)

# def read_externDB(path, db):
#     dbsnp = pd.read_csv(path, sep='\t')
#     print(dbsnp.head())
#     # db.compare_db()



def main():
    # path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/'
    # path_db = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long2.db'
    db = Database(sys.argv[1]) #sys.argv[1]
    # mydb_connection = db.mydb_connection
    # cursor = db.cursor
    read_file(sys.argv[2], db)

    db.close()

main()
