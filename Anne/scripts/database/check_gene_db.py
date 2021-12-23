passwd = ''
# https://www.youtube.com/watch?v=WDEyt2VHpj4
# https://www.youtube.com/watch?v=vR5utJvN4JY

# python3 -m pip install mysql-connector #mysql-connector-python
# python3 -m pip install --upgrade pip

"""
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';

-- -----------------------------------------------------
-- Schema mydb
-- -----------------------------------------------------
-- -----------------------------------------------------
-- Schema internship
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema internship
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS `internship` DEFAULT CHARACTER SET utf8mb4 COLLATE utf8mb4_0900_ai_ci ;
USE `internship` ;

-- -----------------------------------------------------
-- Table `internship`.`project`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `internship`.`project` (
  `ID` INT NOT NULL AUTO_INCREMENT,
  `project_id` VARCHAR(45) NOT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC) VISIBLE,
  UNIQUE INDEX `project_id_UNIQUE` (`project_id` ASC) VISIBLE)
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8mb4
COLLATE = utf8mb4_0900_ai_ci;


-- -----------------------------------------------------
-- Table `internship`.`donor`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `internship`.`donor` (
  `ID` INT NOT NULL AUTO_INCREMENT,
  `donor_ID` VARCHAR(45) NOT NULL,
  `project_ID` INT NOT NULL,
  PRIMARY KEY (`ID`, `project_ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC) VISIBLE,
  UNIQUE INDEX `donor_ID_UNIQUE` (`donor_ID` ASC) VISIBLE,
  INDEX `fk_donor_project_idx` (`project_ID` ASC) VISIBLE,
  CONSTRAINT `fk_donor_project`
    FOREIGN KEY (`project_ID`)
    REFERENCES `internship`.`project` (`ID`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8mb4
COLLATE = utf8mb4_0900_ai_ci;


-- -----------------------------------------------------
-- Table `internship`.`snp`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `internship`.`snp` (
  `ID` INT NOT NULL AUTO_INCREMENT,
  `chr` VARCHAR(45) NULL DEFAULT NULL,
  `pos_start` INT NULL DEFAULT NULL,
  `pos_end` INT NULL DEFAULT NULL,
  `ref` VARCHAR(100) NULL DEFAULT NULL,
  `alt` VARCHAR(100) NULL DEFAULT NULL,
  `genome_version` VARCHAR(45) NULL DEFAULT NULL,
  `depth` INT NULL DEFAULT NULL,
  `GLOC` ENUM('int', 'ext', 'non') NULL DEFAULT NULL,
  `platform` VARCHAR(45) NULL DEFAULT NULL,
  `seq_strategy` VARCHAR(45) NULL DEFAULT NULL,
  `tissue_id` VARCHAR(45) NULL DEFAULT NULL,
  PRIMARY KEY (`ID`),
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC) VISIBLE,
  UNIQUE INDEX `unique_index` (`chr` ASC, `pos_start` ASC, `pos_end` ASC, `ref` ASC, `alt` ASC, `genome_version` ASC, `GLOC` ASC, `platform` ASC, `seq_strategy` ASC) VISIBLE)
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8mb4
COLLATE = utf8mb4_0900_ai_ci;


-- -----------------------------------------------------
-- Table `internship`.`donor_has_snp`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `internship`.`donor_has_snp` (
  `donor_ID` INT NOT NULL,
  `donor_project_ID` INT NOT NULL,
  `snp_ID` INT NOT NULL,
  PRIMARY KEY (`donor_ID`, `donor_project_ID`, `snp_ID`),
  INDEX `fk_donor_has_snp_snp1_idx` (`snp_ID` ASC) VISIBLE,
  INDEX `fk_donor_has_snp_donor1_idx` (`donor_ID` ASC, `donor_project_ID` ASC) VISIBLE,
  CONSTRAINT `fk_donor_has_snp_donor1`
    FOREIGN KEY (`donor_ID` , `donor_project_ID`)
    REFERENCES `internship`.`donor` (`ID` , `project_ID`),
  CONSTRAINT `fk_donor_has_snp_snp1`
    FOREIGN KEY (`snp_ID`)
    REFERENCES `internship`.`snp` (`ID`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8mb4
COLLATE = utf8mb4_0900_ai_ci;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

"""

import mysql.connector
from mysql.connector import connection
import pandas as pd
import glob

def fill_database(mydb_connection, cursor, df):
    """

    :param mydb_connection:
    :param cursor:
    :param df:
    :return:
    """
    # Loop over set of project_ids and add it to the database
    for project_id in list(set(df['project_id'])):
        cursor.execute("""INSERT INTO project (project_id) 
                          VALUES ('%s')""" % (str(project_id)))
        # Get the last ID (private key of the project table) used
        last_id_project = cursor.lastrowid
        print(f"{last_id_project} - {project_id}")
        # Filter dataframe on project_id
        select_project = df.loc[df['project_id'] == project_id]
        print(f"donors: {len(set(select_project['donor_id']))}")
        # Loop over set of donor_ids in (last) project_id and add it to the database
        for donor_id in list(set(select_project['donor_id'])):
            cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                              VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
            # Get the last ID (private key of the donor table) used
            last_id_donor = cursor.lastrowid
            # Filter dataframe on donor_id
            select_donor = select_project[select_project['donor_id'] == donor_id]
            # Loop over rows in dataframe (select_donor)
            for index, row in select_donor.iterrows():
                # See if an SNP already exists with these values
                cursor.execute(
                    """SELECT *
                    FROM snp
                    WHERE chr = '%s' AND pos_start = %s AND pos_end = %s AND 
                    ref = '%s' AND alt = '%s' AND genome_version = '%s' 
                    AND platform = '%s' AND seq_strategy = '%s';""" %
                    (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                     str(row['ref']), str(row['alt']), str(row['genome_version']),
                     str(row['platform']), str(row['seq_strategy'])))
                check_snp = cursor.fetchall()
                # If the SNP does not exist add it to the database
                if not check_snp:
                    cursor.execute("""
                        INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                                     platform, seq_strategy, tissue_id, in_transcript, in_coding, in_exon)
                        VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s', NULL, NULL, NULL)""" %
                                   (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                                    str(row['ref']), str(row['alt']), str(row['genome_version']), int(row['depth']),
                                    str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
                    # Get the last ID (private ket of the snp table) used
                    last_id_snp = cursor.lastrowid
                    # Fill the table donor_has_snp
                    cursor.execute("""
                        INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                        VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(last_id_snp)))
                # If the snp already exists insert the link between the donor and the snp by filling in
                # the donor_has_snp table
                else:
                    # Loop over the snp(s) it corresponds to, if all is well this is always 1 snp.
                    for info in check_snp:
                        # Get ID of the snp
                        id_snp = int(info['ID'])
                        # Check whether the combination donor and snp already exists
                        cursor.execute(
                            """SELECT *
                            FROM donor_has_snp
                            WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s;""" %
                            (int(last_id_project), int(last_id_donor), int(id_snp))
                        )
                        check_donor_snp = cursor.fetchall()
                        # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                        # with this combination
                        if not check_donor_snp:
                            cursor.execute("""
                                INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                                VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(id_snp)))
    # Committing the current transactions
    mydb_connection.commit()



def read_files(path, mydb_connection, cursor):
    """

    :param path:
    :param mydb_connection:
    :param cursor:
    :return:
    """
    path_files = f'{path}*.tsv'
    # Loop over all files
    for fname in glob.glob(path_files):
        print(fname)
        # Read file
        df = pd.read_csv(fname, sep='\t')
        # Drop all duplicates (only depth and tissue_id may differ from all columns)
        df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
        fill_database(mydb_connection, cursor, df)

def check_gene(path_fgene, mydb_connection, cursor):
    gene_df = pd.read_csv(path_fgene, sep='\t')
    print(gene_df.columns)
    for index, row in gene_df.iterrows():
        print('--------------')
        print(index)
        # cursor.execute(
        #             """SELECT *
        #             FROM snp
        #             WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
        #             AND GLOC is NULL;""" %
        #             (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
        # check_snp = cursor.fetchall()
        # for x in check_snp:
        cursor.execute(
                """UPDATE snp
                SET in_transcript = 'true'
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND in_transcript is NULL;""" %
                (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
        mydb_connection.commit()
        cursor.execute(
                """UPDATE snp
                SET in_coding = 'true'
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND in_coding is NULL;""" %
                (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
        mydb_connection.commit()
        exon_start = row['exonStarts'].rstrip(',').split(',')
        exon_end = row['exonEnds'].rstrip(',').split(',')
        # print(exon_start)
        # print(exon_end)
        
        for i in range(int(row['exonCount'])):
            cursor.execute(
                    """UPDATE snp
                    SET in_exon = 'true'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND in_exon is NULL;""" %
                    (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
            mydb_connection.commit()

            #BOOLEAN






def main():
    """

    :return:
    """
    try:
        # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
        mydb_connection = mysql.connector.connect(host='127.0.0.1', user='root', passwd=passwd, database='internship')
        # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
        cursor = mydb_connection.cursor(dictionary=True)
        path_fgene = 'D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed'
        path = 'D:/Hanze_Groningen/STAGE/db/files/'
        
        # read_files(path, mydb_connection, cursor)
        check_gene(path_fgene, mydb_connection, cursor)

        # cursor.execute('SELECT * FROM donor_has_snp')
        # for x in cursor:
        #     print(f"{x['donor_project_ID']} - {x['donor_ID']} - {x['snp_ID']}")
    except mysql.connector.Error as e:
        print("Error while connecting to sqlite", e)
    finally:
        if mydb_connection:
            mydb_connection.close()
            print("The SQLite connection is closed")


main()
print('END')