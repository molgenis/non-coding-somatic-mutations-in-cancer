#!/usr/bin/env python3

# Imports
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def set_dosages(db):
    """
    Add a new value to the database: dosages.
    Set dosages to a new value.
    :param db:  The database object
    :return:
    """
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `dosages` FLOAT NULL DEFAULT NULL
                    """)
    # set dosages to the correct value
    db.cursor.execute(
        """UPDATE donor_has_snp 
                SET dosages = (CAST(mutant_allele_read_count AS REAL) / CAST(total_read_count AS REAL))
                WHERE total_read_count > 0;""")
    # Add to database
    db.mydb_connection.commit()


def set_GT(db, table):
    """
    Add and set GT (genotype) to homozygous alt (2), heterozygous (1) or homozygous ref (0).
    GT remains Nan (or none) when there are no dosages.
    :param db:    The database object
    :param table: Name of the table
    :return:
    """
    # add dosages
    db.cursor.execute("""
                    ALTER TABLE %s
                    ADD `GT` INT NULL DEFAULT NULL
                    """ % table)
    # set GT to the correct value
    # Homozygoot alt
    db.cursor.execute(
        """UPDATE %s 
                SET GT = 2
                WHERE dosages >= 0.8;""" % table)
    # Heterozygoot 
    db.cursor.execute(
        """UPDATE %s 
                SET GT = 1
                WHERE dosages < 0.8 AND dosages > 0.2;""" % table)
    # Homozygoot ref
    db.cursor.execute(
        """UPDATE %s 
                SET GT = 0
                WHERE dosages <= 0.2;""" % table)
    # Add to database
    db.mydb_connection.commit()


def set_GT2(db, table):
    """
    Add and set GT2 (genotype).
    The Nan (or none) values are all set to 0 (homozygous ref).
    :param db:    The database object
    :param table: Name of the table
    :return:
    """
    # add dosages
    db.cursor.execute("""
                    ALTER TABLE %s
                    ADD `GT2` INT NULL DEFAULT NULL
                    """ % table)
    # set GT2 to the correct value
    # Make all values of GT2 equal to GT
    db.cursor.execute(
        """UPDATE %s 
                SET GT2 = GT;""" % table)
    # Set all Nan (or none) values to 0 (homozygous ref)
    db.cursor.execute(
        """UPDATE %s 
                SET GT2 = 0
                WHERE GT IS NULL;""" % table)
    db.mydb_connection.commit()


def sum_dosage_GT(db):
    """
    Add new table based on 'donor_had_snp'.
    :param db:  The database object
    :return:
    """
    # Create new table
    db.cursor.execute("""
        CREATE TABLE sum_dosage_GT AS
            SELECT donor_has_snp.snp_ID, donor_has_snp.donor_ID, donor_has_snp.donor_project_ID,
                    (CAST(SUM(donor_has_snp.mutant_allele_read_count) AS REAL) / 
                    CAST(SUM(donor_has_snp.total_read_count) AS REAL)) AS "dosages",
                    SUM(donor_has_snp.total_read_count) AS "total_read_count_sum", 
                    SUM(donor_has_snp.mutant_allele_read_count) AS "mutant_allele_read_count_sum",
                    COUNT(donor_has_snp.total_read_count) AS "number_snps"
            FROM donor_has_snp
            WHERE donor_has_snp.total_read_count > 0
            GROUP BY donor_has_snp.donor_ID, donor_has_snp.snp_ID;
    """)
    # https://database.guide/add-a-foreign-key-to-an-existing-table-in-sqlite/
    db.cursor.execute("""
        PRAGMA foreign_keys = OFF;
    """)
    db.cursor.execute("""
        BEGIN TRANSACTION;
    """)
    # Create new table with forend keys
    db.cursor.execute("""
        CREATE TABLE sum_dosage_GT_new( 
                    `snp_ID` INT NOT NULL,
                    `donor_ID` INT NOT NULL,
                    `donor_project_ID` INT NOT NULL,
                    `dosages` INT NULL DEFAULT NULL,
                    `total_read_count_sum` INT NULL DEFAULT NULL,
                    `mutant_allele_read_count_sum` INT NULL DEFAULT NULL,
                    `number_snps` INT NULL DEFAULT NULL,
                    
                    CONSTRAINT `fk_sum_dosage_GT_new`
                        FOREIGN KEY (`donor_ID`, `donor_project_ID`, `snp_ID`)
                        REFERENCES `donor_has_snp` (`donor_ID`, `donor_project_ID`, `snp_ID`)
                );
    """)
    db.cursor.execute("""
        INSERT INTO sum_dosage_GT_new SELECT * FROM sum_dosage_GT;
    """)
    db.cursor.execute("""
        DROP TABLE sum_dosage_GT;
    """)
    db.cursor.execute("""
        ALTER TABLE sum_dosage_GT_new RENAME TO sum_dosage_GT;
    """)
    db.cursor.execute("""
        COMMIT;
    """)
    db.cursor.execute("""
        PRAGMA foreign_keys = ON;
    """)
    # Committing the current transactions
    db.mydb_connection.commit()


def main():
    # Call get_config
    config = get_config('gearshift')
    # Path of the database
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    set_dosages(db)
    set_GT(db, 'donor_has_snp')
    set_GT2(db, 'donor_has_snp')
    sum_dosage_GT(db)
    set_GT(db, 'sum_dosage_GT')
    set_GT2(db, 'sum_dosage_GT')


if __name__ == '__main__':
    main()
