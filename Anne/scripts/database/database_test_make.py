import sqlite3
import glob
import pandas as pd
import sys


def create_db(cursor):
    """

    :param cursor:
    :return:
    """
    # Create tables
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `project`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `project`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `project_ID` VARCHAR(45) NOT NULL  
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `donor`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `donor_ID` VARCHAR(45) NOT NULL,
        `project_ID` INT NOT NULL,
        CONSTRAINT `fk_donor_project`
            FOREIGN KEY (`project_ID`)
            REFERENCES `project` (`ID`)
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `snp`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `chr` VARCHAR(45) NULL DEFAULT NULL,
        `pos_start` INT NULL DEFAULT NULL,
        `pos_end` INT NULL DEFAULT NULL,
        `ref` VARCHAR(1000) NULL DEFAULT NULL,
        `alt` VARCHAR(1000) NULL DEFAULT NULL,
        `genome_version` VARCHAR(45) NULL DEFAULT NULL,
        `depth` INT NULL DEFAULT NULL,
        `in_transcript` BOOLEAN DEFAULT(FALSE),
        `in_coding` BOOLEAN DEFAULT(FALSE),
        `in_exon` BOOLEAN DEFAULT(FALSE),
        `platform` VARCHAR(45) NULL DEFAULT NULL,
        `seq_strategy` VARCHAR(45) NULL DEFAULT NULL,
        `tissue_id` VARCHAR(45) NULL DEFAULT NULL,
        UNIQUE (`chr`, `pos_start`, `pos_end`, `ref`, `alt`, `genome_version`, `platform`, `seq_strategy`) 
        ON CONFLICT REPLACE
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `donor_has_snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor_has_snp`(
        `donor_ID` INT NOT NULL,
        `donor_project_ID` INT NOT NULL,
        `snp_ID` INT NOT NULL,
        CONSTRAINT `fk_donor_has_snp_donor1`
            FOREIGN KEY (`donor_ID`, `donor_project_ID`)
            REFERENCES `donor` (`ID`, `project_ID`)
        CONSTRAINT `fk_donor_has_snp_snp1`
            FOREIGN KEY (`snp_ID`)
            REFERENCES `snp` (`ID`)
    )
    """)
    return cursor


def main():
    """

    :return:
    """
    path_db = f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/{sys.argv[1]}.db'
    try:
        # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
        mydb_connection = sqlite3.connect(path_db)
        # Setting row_factory property of connection object to
        # sqlite3.Row(sqlite3.Row is an implementation of row_factory)
        mydb_connection.row_factory = sqlite3.Row
        # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
        cursor = mydb_connection.cursor()
        cursor = create_db(cursor)
        
    except sqlite3.Error as er:
        print("Error while connecting to sqlite", er)
    finally:
        if mydb_connection:
            mydb_connection.close()
            print("The SQLite connection is closed")


main()
print(f'END {sys.argv[1]}')
