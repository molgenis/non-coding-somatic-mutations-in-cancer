import sqlite3
import glob
from sqlite3.dbapi2 import connect
import pandas as pd
import numpy as np


def create_db(name):
    # returns a connection object that we will use to interact with the SQLite database held in the file test.db
    connection = sqlite3.connect(f'{name}.db')
    print(connection.total_changes)
    # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
    cursor = connection.cursor()
    # Create tables
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `internship`.`project`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `project`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `project_ID` VARCHAR(45) NOT NULL  
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `internship`.`donor`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor`(
        `ID` INT NOT NULL PRIMARY KEY UNIQUE,
        `donor_ID` VARCHAR(45) NOT NULL,
        `project_ID` INT NOT NULL,
        CONSTRAINT `fk_donor_project`
            FOREIGN KEY (`project_ID`)
            REFERENCES `project` (`ID`)
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `internship`.`snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `snp`(
        `ID` INT NOT NULL PRIMARY KEY UNIQUE,
        `chr` VARCHAR(45) NULL DEFAULT NULL,
        `pos_start` INT NULL DEFAULT NULL,
        `pos_end` INT NULL DEFAULT NULL,
        `ref` VARCHAR(1000) NULL DEFAULT NULL,
        `alt` VARCHAR(1000) NULL DEFAULT NULL,
        `genome_version` VARCHAR(45) NULL DEFAULT NULL,
        `depth` INT NULL DEFAULT NULL,
        `GLOC` TEXT CHECK( `GLOC` IN ('int','ext','non') ) NULL DEFAULT NULL,
        `platform` VARCHAR(45) NULL DEFAULT NULL,
        `seq_strategy` VARCHAR(45) NULL DEFAULT NULL,
        `tissue_id` VARCHAR(45) NULL DEFAULT NULL,
        UNIQUE (`chr`, `pos_start`, `pos_end`, `ref`, `alt`, `genome_version`, `GLOC`, `platform`, `seq_strategy`) ON CONFLICT REPLACE
    )
    """)
    cursor.execute("""
    -- -----------------------------------------------------
    -- Table `internship`.`donor_has_snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor_had_snp`(
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
    return connection, cursor

def fill_database(connection, cursor, df):
    for project_id in list(set(df['project_id'])):
        cursor.execute("""INSERT INTO project (project_id) VALUES ('%s')""" % (str(project_id)))
    connection.commit()

def read_files(path, connection, cursor):
    path_files = f'{path}*.tsv'
    # Loop over all files
    for fname in glob.glob(path_files):
        print(fname)
        # Read file
        df = pd.read_csv(fname, sep='\t')
        # Drop all duplicates (die hetzelfde zijn behlave de depth deze mag verschillen)
        df = df.drop_duplicates(subset=df.columns.difference(['depth']))
        fill_database(connection, cursor, df)
        



connection, cursor = create_db('test')
read_files('D:/Hanze_Groningen/STAGE/db/files/', connection, cursor)

cursor.execute('SELECT * FROM project')
for x in cursor:
    print(x)
