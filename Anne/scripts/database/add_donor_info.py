#!/usr/bin/env python3
from pickle import FALSE, TRUE
import pandas as pd
import sys
import numpy as np
from Database import Database


def add_columns(cursor, mydb_connection):
    # sex --> false = man en later passen we zit aan als het een vrouw is naar 1=true
    # vital_status --> false = alive en later passen we dit aan wanneer nodig 1=true=deceased
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `sex` BOOLEAN DEFAULT(FALSE)
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `vital_status` BOOLEAN DEFAULT(FALSE)
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `age_at_diagnosis` INT NULL DEFAULT NULL
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `age_at_enrollment` INT NULL DEFAULT NULL
                    """)
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `age_at_last_followup` INT NULL DEFAULT NULL
                    """)
    # cursor.execute(f"""
    #                 ALTER TABLE donor
    #                 ADD `relapse_interval` INT NULL DEFAULT NULL
    #                 """)
    cursor.execute(f"""
                    ALTER TABLE donor
                    ADD `survival_time` INT NULL DEFAULT NULL
                    """)
    mydb_connection.commit()

def drop_col(cursor, mydb_connection):
    cursor.execute(f"""
                    ALTER TABLE donor DROP COLUMN `relapse_interval`;
                    """)
    mydb_connection.commit()
    

def add_info(cursor, mydb_connection, donor_df):
    # Loop over rows in eQTL file
    for index, row in donor_df.iterrows():
        if row['donor_sex'] == 'male':
            sex = str('FALSE')
        elif row['donor_sex'] == 'female':
            sex = str('TRUE')

        if row['donor_vital_status'] == 'deceased':
            vital_status = str('FALSE')
        elif row['donor_vital_status'] == 'alive':
            vital_status = str('TRUE')
        print('-----')
        print(sex)
        print(vital_status)
        print(row['donor_age_at_diagnosis'])
        print(row['donor_age_at_enrollment'])
        print(row['donor_age_at_last_followup'])
        # print(row['donor_relapse_interval'])
        print(row['donor_survival_time'])
        print(row['icgc_donor_id'])
        
        # Add ID_eQTL and eQTL to snp that matches the matches in chr, pos_start, pos_end, ref and alt.
        #, project_ID = %s, relapse_interval = %s, 
        cursor.execute(
                """UPDATE `donor`
                    SET sex = %s, vital_status = %s, age_at_diagnosis = %s,
                    age_at_enrollment = %s, age_at_last_followup = %s,
                    survival_time = %s
                    WHERE donor_ID = '%s';""" %
                (sex, vital_status, int(row['donor_age_at_diagnosis']), 
                int(row['donor_age_at_enrollment']), int(row['donor_age_at_last_followup']),
                int(row['donor_survival_time']), str(row['icgc_donor_id'])))
    # Add to database
    mydb_connection.commit()

def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db"
    # Path of the genes and there positions
    donor_path = "D:/Hanze_Groningen/STAGE/metadata CANCER/ALL-US/donor.tsv"
    # Read gene file
    donor_df = pd.read_csv(donor_path, sep='\t')  # sys.argv[2], sep='\t')
    print(len(donor_df))
    # Database connection
    db = Database(path_db)  # sys.argv[1]
    # Call ...
    # add_columns(db.cursor, db.mydb_connection)
    drop_col(db.cursor, db.mydb_connection)
    add_info(db.cursor, db.mydb_connection, donor_df)
    # Close connections
    db.close()


if __name__ == '__main__':
    main()





