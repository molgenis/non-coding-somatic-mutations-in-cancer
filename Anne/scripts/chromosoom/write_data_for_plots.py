#!/usr/bin/env python3
import sys
import os

# Also takes the folder 1 higher, so that I can do the import after
# sys.path.append("..")
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

def IN_OUT(cursor, info, column_name, file, true_false):
    cursor.execute("""
            SELECT *
            FROM snp
            WHERE ID = %s AND seq_strategy = 'WGS' AND %s = %s
            """%
            (int(info['snp_ID']), str(column_name), int(true_false)))
    snp_info = cursor.fetchall()
    # Loop over snps
    for inf_snp in snp_info:
        # Add the information in the way that R needs to plot later.
        file.write(f"{inf_snp['ID']}\tchr{inf_snp['chr']}\t{inf_snp['pos_start']}\t{inf_snp['pos_end']}\n")
    
    return file



def make_karyoploteR(db, path_files, column_name):
    """
    Create the files that karyoploteR uses to make pictures
    :param db:              the database object
    :param path_files:      the path to save the files
    :return:
    """
    # Find all donors
    db.cursor.execute("""
                    SELECT *
                    FROM donor
                    """)
    all_donors = db.cursor.fetchall()
    # Loop over all donors
    for don in all_donors:
        # Find all snps IDz where donor_ID is equal to donor and where donor_project_ID is equal to project ID.
        db.cursor.execute(f"""
                        SELECT *
                        FROM donor_has_snp
                        WHERE donor_ID = {don['ID']} AND donor_project_ID = {don['project_ID']}
                        """)
        all_snps = db.cursor.fetchall()
        # Open and create file
        f = open(f"{path_files}{don['project_ID']}_{don['donor_ID']}.bed", "a")
        f_IN = open(f"{path_files}{don['project_ID']}_{don['donor_ID']}_IN_{column_name}.bed", "a")
        f_OUT = open(f"{path_files}{don['project_ID']}_{don['donor_ID']}_OUT_{column_name}.bed", "a")
        # Loop over all snps
        for info in all_snps:
            # Find all snps that match the SNP id and have used WGS
            db.cursor.execute("""
                    SELECT *
                    FROM snp
                    WHERE ID = %s AND seq_strategy = 'WGS'
                    """%
                    (info['snp_ID']))
            snp_info = db.cursor.fetchall()
            # Loop over snps
            for inf_snp in snp_info:
                # Add the information in the way that R needs to plot later.
                f.write(f"{inf_snp['ID']}\tchr{inf_snp['chr']}\t{inf_snp['pos_start']}\t{inf_snp['pos_end']}\n")

            f_IN = IN_OUT(db.cursor, info, column_name, f_IN, 1)
            f_OUT = IN_OUT(db.cursor, info, column_name, f_OUT, 0)            

        f.close()
        f_IN.close()
        f_OUT.close()


def main():
    """

    """
    config = get_config('gearshift')
    # Path to database
    path_db = "D:/Hanze_Groningen/STAGE/TEST_DEL/test3.db"
    # Path to where the plots should be saved
    path_to = 'D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/'
    path_files = f'{path_to}other/'
    # Makes folder if it doesn't exists
    if not os.path.exists(path_files):
        os.makedirs(path_files)
    # Makes Database object
    db = Database(path_db)
    # Calls make_karyoploteR
    make_karyoploteR(db, path_files, 'in_transcript')
    db.close()


if __name__ == '__main__':
    main()
