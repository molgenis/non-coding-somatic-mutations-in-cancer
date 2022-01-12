import sqlite3
import glob
# import pandas as pd
import sys
import io
import os

from Database import Database



def make_karyoploteR(db, path_files):
    db.cursor.execute("""
                    SELECT *
                    FROM donor
                    """)
    all_donors = db.cursor.fetchall()
    for don in all_donors:
        db.cursor.execute(f"""
                        SELECT *
                        FROM donor_has_snp
                        WHERE donor_ID = {don['ID']} AND donor_project_ID = {don['project_ID']}
                        """)
        all_snps = db.cursor.fetchall()
        f = open(f"{path_files}{don['project_ID']}_{don['donor_ID']}.bed", "a")
        for info in all_snps:
            db.cursor.execute(f"""
                    SELECT *
                    FROM snp
                    WHERE ID = {info['snp_ID']} AND seq_strategy = 'WGS'
                    """)
            snp_info = db.cursor.fetchall()                                
            for inf_snp in snp_info:
                f.write(f"{inf_snp['ID']}\tchr{inf_snp['chr']}\t{inf_snp['pos_start']}\t{inf_snp['pos_end']}\n")
        f.close()



def main():
    path_db = "D:/Hanze_Groningen/STAGE/TEST_DEL/test3.db"
    path_files = 'D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/'
    db = Database(path_db)
    make_karyoploteR(db, path_files)


if __name__ == '__main__':
    main()