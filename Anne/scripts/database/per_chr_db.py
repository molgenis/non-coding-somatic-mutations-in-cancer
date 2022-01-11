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



def fill_project_donor(df, db):
    # db = Database(path_db) #sys.argv[1]
    mydb_connection = db.mydb_connection
    cursor = db.cursor
    # Loop over set of project_ids and add it to the database
    for project_id in sorted(list(set(df['project_id'])), key=str.lower):
        print(project_id)
        cursor.execute("""INSERT INTO project (project_ID) 
                        VALUES ('%s')""" % (str(project_id)))
        # Committing the current transactions
        mydb_connection.commit()
        # Get the last ID (private key of the project table) used
        last_id_project = cursor.lastrowid
        print(f"{last_id_project} - {project_id}")
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        print(f"donors: {len(set(select_project['donor_id']))}")
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in sorted(list(set(select_project['donor_id'])), key=str.lower):
            cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                            VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
            # Committing the current transactions
            mydb_connection.commit()
            # # Get the last ID (private key of the donor table) used
            # last_id_donor = cursor.lastrowid
            # # Filter dataframe on donor_id
            # select_donor = select_project[select_project['donor_id'] == donor_id]

            # print(mp.cpu_count())
            # df_shuffled = select_donor.sample(frac=1)
            # df_splits = np.array_split(df_shuffled, mp.cpu_count())
            # arg_multi_list = []
            # for df_s in df_splits:
            #     arg_multi_list.append((df_s, path_db, last_id_project, last_id_donor))

            # pool = Pool(processes=mp.cpu_count())
            # pool.starmap(func=fill_snp, iterable=arg_multi_list)
            # pool.close()
            # pool.join()

    db.close()




def read_file(path, db):
    """
    
    :param path: 
    :param db: 
    :return: 
    """
    df = pd.read_csv(path, sep='\t')
    # Drop all duplicates (only depth and tissue_id may differ from all columns)
    df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
    # df ,sys.argv[1]
    fill_project_donor(df, db)


def main():
    # path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/data_db/'
    # path_db = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long2.db'
    
    db = Database(sys.argv[1]) #sys.argv[1]
    # mydb_connection = db.mydb_connection
    # cursor = db.cursor
    read_file(sys.argv[2], db)

    db.close()

if __name__ == '__main__':
    main()