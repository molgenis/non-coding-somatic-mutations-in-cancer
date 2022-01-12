#!/usr/bin/env python3
import sys
import os


# Also takes the folder 1 higher, so that I can do the import after
sys.path.append("..")
from database.Database import Database



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
    path_to = 'D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/'
    path_files = f'{path_to}other/'
    if not os.path.exists(path_files):
        os.makedirs(path_files)
    db = Database(path_db)
    make_karyoploteR(db, path_files)


if __name__ == '__main__':
    main()