#!/usr/bin/env python3

# Imports
import sys
sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def get_min_max_snpID(db):
    """
    Finds the highest and lowest snp ID.
    :param db:  The database object
    :return: max_snp_id: The highest snp id
             min_snp_id: The lowest snp id
    """
    # MAX
    db.cursor.execute("""
                    SELECT MAX(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        max_snp_id = res[0]
    # MIN
    db.cursor.execute("""
                    SELECT MIN(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        min_snp_id = res[0]
    return max_snp_id, min_snp_id


def add_value(db):
    """
    Adds values (AF, AF_breast and AF_nonbreast) to the database (table snp).
    :param db:  The database object
    :return:
    """
    # Add AF
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `AF` INT NULL DEFAULT NULL
                    """)
    # Add AF_breast
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `AF_breast` INT NULL DEFAULT NULL
                    """)
    # Add AF_nonbreast
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `AF_nonbreast` INT NULL DEFAULT NULL
                    """)
    # Committing the current transactions
    db.mydb_connection.commit()


def cal_AF(db):
    """
    Calculate AF
    :param db:  The database object
    :return:
    """
    db.cursor.execute(
        """UPDATE snp
            SET AF = cal_AF.AF
            FROM (SELECT (CAST(SUM(sum_dosage_GT.GT2) AS REAL) / (CAST(COUNT(sum_dosage_GT.donor_ID) AS REAL) * 2)) 
                    AS AF, 
                    FROM sum_dosage_GT, snp
                    WHERE sum_dosage_GT.snp_ID = snp.ID AND (sum_dosage_GT.GT2 = 0 OR sum_dosage_GT.GT2 = 1 
                    OR sum_dosage_GT.GT2 = 2)
                    GROUP BY sum_dosage_GT.snp_ID
                    ORDER BY sum_dosage_GT.snp_ID) AS cal_AF;""")
    # Add to database
    db.mydb_connection.commit()


def main():
    # Call get_config
    config = get_config('gearshift')
    # Path of the database
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    # Call add_value
    add_value(db)
    # Call cal_AF
    cal_AF(db)


if __name__ == '__main__':
    main()
