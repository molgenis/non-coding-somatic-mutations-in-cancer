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
  UNIQUE INDEX `ID_UNIQUE` (`ID` ASC) VISIBLE)
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
import pandas as pd

mydb = mysql.connector.connect(host='127.0.0.1', user='root', passwd=passwd, database='internship')
mycursor = mydb.cursor()

mycursor.execute("SHOW TABLES")
for x in mycursor:
    print(x)

path = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/db/BOCA-UK_db.tsv"
df = pd.read_csv(path, sep='\t')
print(df.columns)

def per_row():
  for project_id in list(set(df['project_id'])):
      mycursor.execute("""INSERT INTO project (project_id)
                    VALUES ('%s')""" % (str(project_id)))
      last_id_project = mycursor.lastrowid
      print(last_id_project)  
      select_project = df.loc[df['project_id'] == project_id]
      print(len(set(select_project['donor_id'])))  
      for donor_id in list(set(select_project['donor_id']))[:2]:        
          mycursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                    VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
          last_id_donor = mycursor.lastrowid
          

          select_donor = df.loc[df['donor_id'] == donor_id]
          for index, row in df.iterrows():
            if index < 20:
              print(index)
              # TODO KIJKEN HOE JE DEZELFDE SNPS MOET OPSLAAN
              mycursor.execute("""INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, GLOC, platform, seq_strategy, tissue_id)
                    VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s', '%s')""" % (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                    str(row['ref']), str(row['alt']), str(row['genome_version']), int(row['depth']), str('int'), str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
              last_id_snp = mycursor.lastrowid
              mycursor.execute("""INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                    VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(last_id_snp)))
      mydb.commit()
      mycursor.execute('SELECT * FROM donor')
      for x in mycursor:
          print(x)

def per_tabel():
  print()

per_row()


print('END')