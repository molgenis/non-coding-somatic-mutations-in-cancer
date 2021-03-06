#!/usr/bin/env python3
import pandas as pd
import sys

from Database import Database


def create_db(db):
    """
    Create the database
    :param db:  the database object
    :return:
    """
    # Create tables
    db.cursor.execute("""
    -- -----------------------------------------------------
    -- Table `project`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `project`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `project_ID` VARCHAR(45) NOT NULL  
    )
    """)
    db.cursor.execute("""
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
    db.cursor.execute("""
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
    db.cursor.execute("""
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


def main():
    """

    :return:
    """
    # Make Database object
    db = Database(sys.argv[1])
    # Call create_db
    create_db(db)
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
