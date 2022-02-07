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
        `project_ID` VARCHAR(45) NOT NULL,
        'cancer' VARCHAR(45) NULL DEFAULT NULL
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
        `sex` BOOLEAN DEFAULT(FALSE),
        `vital_status` BOOLEAN DEFAULT(FALSE), 
        `age_at_diagnosis` INT NULL DEFAULT NULL, 
        `age_at_last_followup` INT NULL DEFAULT NULL,
        `disease_status_last_followup` VARCHAR(45) DEFAULT NULL,
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
        `platform` VARCHAR(45) NULL DEFAULT NULL,
        `seq_strategy` VARCHAR(45) NULL DEFAULT NULL,
        UNIQUE (`chr`, `pos_start`, `pos_end`, `ref`, `alt`, `genome_version`, `platform`, `seq_strategy`) 
        ON CONFLICT REPLACE
    )
    """)
    print('HAI')
    db.cursor.execute("""
    -- -----------------------------------------------------
    -- Table `tissue`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `tissue`(
        `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
        `specimen_type` VARCHAR(100) NOT NULL,
        `type` BOOLEAN DEFAULT(FALSE)
    )
    """)
    print('HALLO')
    db.cursor.execute("""
    -- -----------------------------------------------------
    -- Table `donor_has_snp`
    -- -----------------------------------------------------
    CREATE TABLE IF NOT EXISTS `donor_has_snp`(
        `donor_ID` INT NOT NULL,
        `donor_project_ID` INT NOT NULL,
        `snp_ID` INT NOT NULL,
        `tissue_ID`  INT NOT NULL,
        `specimen_id` VARCHAR(45) NOT NULL,
        CONSTRAINT `fk_donor_has_snp_donor1`
            FOREIGN KEY (`donor_ID`, `donor_project_ID`)
            REFERENCES `donor` (`ID`, `project_ID`)
        CONSTRAINT `fk_donor_has_snp_snp1`
            FOREIGN KEY (`snp_ID`)
            REFERENCES `snp` (`ID`)
        CONSTRAINT `fk_donor_has_snp_tissue1`
            FOREIGN KEY (`tissue_ID`)
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
