from Database import Database
import sys


def set_dosages(db):
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `dosages` FLOAT NULL DEFAULT NULL
                    """)
    print('hoi')
    # set dosages to the correct value
    db.cursor.execute(
            """UPDATE donor_has_snp 
                SET dosages = (CAST(mutant_allele_read_count AS REAL) / CAST(total_read_count AS REAL));""")
    # Add to database
    db.mydb_connection.commit()
    
    # Check if dosages really changed
    # db.cursor.execute("""
    #                 SELECT *
    #                 FROM 'donor_has_snp'
    #                 """ )
    # check = db.cursor.fetchall()
    # print(len(check))
    # print('CHECK')
    # for ch in check:
    #     if not isinstance(ch['mutant_allele_read_count'], type(None)):
    #         print(f"{ch['mutant_allele_read_count']} / {ch['total_read_count']} = {ch['dosages']}")

def set_GT(db):
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `GT` INT NULL DEFAULT NULL
                    """)
    print('add done')
    # set GT to the correct value
    # Homozygoot alt
    db.cursor.execute(
            """UPDATE donor_has_snp 
                SET GT = 2
                WHERE dosages >= 0.9;""")
    print('homozyg alt done')
    # Heterozygoot 
    db.cursor.execute(
            """UPDATE donor_has_snp 
                SET GT = 1
                WHERE dosages < 0.9 AND dosages > 0.1;""")
    print('hetero done')
    # Homozygoot ref
    db.cursor.execute(
            """UPDATE donor_has_snp 
                SET GT = 0
                WHERE dosages <= 0.1;""")
    print('homozyg ref done')
    # Add to database
    db.mydb_connection.commit()
    
    # Check if dosages really changed
    db.cursor.execute("""
                    SELECT *
                    FROM 'donor_has_snp'
                    """ )
    check = db.cursor.fetchall()
    print(len(check))
    print('CHECK')
    for index, ch in enumerate(check):
        if index <= 10:
            if not isinstance(ch['mutant_allele_read_count'], type(None)):
                print(f"{ch['mutant_allele_read_count']} / {ch['total_read_count']} = {ch['dosages']} --> {ch['GT']}")
    


def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    set_dosages(db)
    set_GT(db)
      



if __name__ == '__main__':
    main()