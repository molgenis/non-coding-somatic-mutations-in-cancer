#!/usr/bin/env python3
import pandas as pd
import sys
import numpy as np
from Database import Database


def add_columns(cursor, mydb_connection):
    cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `tissue` BOOLEAN DEFAULT(FALSE)
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `tissu_type` VARCHAR(100) NULL DEFAULT NULL
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `icgc_specimen_id` VARCHAR(45) NULL DEFAULT NULL
                    """)
    mydb_connection.commit()


def add_info(cursor, mydb_connection, hc_tum_df):
    cursor.execute("""SELECT *
                    FROM project""")
    project = cursor.fetchall()
    for pro in project:
        # print('-------')
        # print(pro['ID'])
        cursor.execute("""SELECT *
                          FROM donor
                          WHERE project_ID = %s"""% 
                          (int(pro['ID'])))
        donors = cursor.fetchall()
        select_project = hc_tum_df[hc_tum_df['project_code'] == pro['ID']]
        for don in donors:
            cursor.execute("""SELECT *
                            FROM donor_has_snp
                            WHERE donor_project_ID = %s AND donor_ID = %s"""% 
                            (int(pro['ID']), int(don['ID'])))
            all_snps = cursor.fetchall()
            select_donor = select_project[select_project['icgc_donor_id'] == don['ID']]
            for snp in all_snps:
                cursor.execute("""SELECT *
                                FROM snp
                                WHERE ID = %s AND seq_strategy = 'WGS'"""% 
                                (int(snp['snp_ID'])))
                info_snps = cursor.fetchall()
                print('###############')
                for info in info_snps:
                    print('-------')
                    print(f"{info['chr']} - {info['pos_start']} - {info['pos_end']} - {info['ref']} - {info['alt']}")

                    # TODO
                    # # FILTER FILE MET CHR POS REF ALT
                    # # GET icgc_specimen_id
                    # # Zoek deze icgc_specimen_id weer op in all_specimen.tsv
                    # select_specimen_id = TEST[TEST['icgc_specimen_id'] == OTHER]
                    # # Get specimen_type
                    # specimen_type = select_specimen_id['specimen_type']
                    # kern = specimen_type.split('-')[0].strip()
                    # if kern == 'Normal':
                    #     # Alter values in donor_has_snp
                    #     # icgc_specimen_id kan ook als extra kolom in snp
                    #     cursor.execute("""UPDATE `donor_has_snp`
                    #                 SET tissue = %s, tissu_type = %s, icgc_specimen_id = %s
                    #                 WHERE donor_ID = '%s' AND donor_project_ID = %s AND snp_ID = %s;""" %
                    #             (TISSUE, TISSUE_TYPE, SPECIMEN_ID, 
                    #             int(don['ID']), int(pro['ID']), int(info['ID'])))
                    # else:
                    #     # Alter values in donor_has_snp
                    #     # icgc_specimen_id kan ook als extra kolom in snp
                    #     cursor.execute("""UPDATE `donor_has_snp`
                    #                 SET tissue = %s, tissu_type = %s, icgc_specimen_id = %s
                    #                 WHERE donor_ID = '%s' AND donor_project_ID = %s AND snp_ID = %s;""" %
                    #             (TISSUE, TISSUE_TYPE, SPECIMEN_ID, 
                    #             int(don['ID']), int(pro['ID']), int(info['ID'])))

        

    # for specimen_type in list(set(hc_tum_df['specimen_type'])):
    #     select_project = hc_tum_df.loc[hc_tum_df['specimen_type'] == specimen_type]
    #     # print('---')
    #     print(specimen_type)
    #     # print(len(list(set(select_project['icgc_donor_id']))))
    #     # for donor_id in list(set(select_project['icgc_donor_id'])):
    #     #     print(donor_id)



def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db"
    # Path of the genes and there positions
    hc_tum_path = "D:/Hanze_Groningen/STAGE/metadata CANCER/WGS/all_specimen.tsv" #"E:/STAGE/WGS/all_specimen.tsv"
    # Read gene file
    hc_tum_df = pd.read_csv(hc_tum_path, sep='\t')  # sys.argv[2], sep='\t')
    print(len(hc_tum_df))
    # Database connection
    db = Database(path_db)  # sys.argv[1]
    # Call ...
    # add_columns(db.cursor, db.mydb_connection)
    # drop_col(db.cursor, db.mydb_connection)
    add_info(db.cursor, db.mydb_connection, hc_tum_df)
    # Close connections
    db.close()


if __name__ == '__main__':
    main()





