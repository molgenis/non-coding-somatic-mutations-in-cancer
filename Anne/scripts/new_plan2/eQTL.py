#!/usr/bin/env python3
import pandas as pd
import sys
from Database import Database


def add_value(db, ID_eQT, eQT, close_eQT):
    """
    Adds values (ID_eQT, eQT, close_eQT) to the database (table snp).
    :param db:  The database object
    :param ID_eQT:  String. This is the name of the parameter from the database that is set to an ID 
                    when that snp exactly matches an eQT.
    :param eQT:  String. This is the name of the parameter in the database that will 
                 be set to true when a snp matches an eQT exactly
    :param close_eQT:  String. This is the name of the parameter in the database that is set 
                       to true when a snp falls within a certain region (before or after) a known eQT
    :return:
    """
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` VARCHAR(45) NULL DEFAULT NULL
                    """ %
                    (ID_eQT))
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` BOOLEAN DEFAULT(FALSE)
                    """ %
                    (eQT))
    db.cursor.execute("""
                    ALTER TABLE snp
                    ADD `%s` BOOLEAN DEFAULT(FALSE)
                    """ %
                    (close_eQT))

def set_value(db, ID_eQT, eQT, row):
    """
    The parameters ID_eQT and eQT are set to either an ID or true when a snp exactly matches an eQT in position, ref and alt
    :param db:  The database object
    :param ID_eQT:  String. This is the name of the parameter from the database that is set to an ID 
                    when that snp exactly matches an eQT.
    :param eQT:  String. This is the name of the parameter in the database that will 
                 be set to true when a snp matches an eQT exactly
    :param row:  Line from the eQT file with information about 1 eQT
    :return:
    """ #TODO also alt and ref????
    # Add ID_eQTL and eQTL to snp that matches the matches in chr, pos_start, pos_end, ref and alt.
    db.cursor.execute(
            """UPDATE snp
                SET %s = '%s', %s = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                AND ref = '%s' AND alt = '%s';""" %
            (ID_eQT, str(row['SNP']), eQT, str(row['SNPChr']), int(row['SNPPos']), 
            int(row['SNPPos']), str(row['AssessedAllele']), str(row['OtherAllele'])))


def set_value_close_to(db, row, close_eQT, region):
    """
    The parameter close_eQT is set to true when a snp falls within the region (before or after) of an eQT
    :param db:  The database object
    :param row:  Line from the eQT file with information about 1 eQT
    :param close_eQT:  String. This is the name of the parameter in the database that is set 
                       to true when a snp falls within a certain region (before or after) a known eQT
    :param region:  Region (front or back) an eQT, within which snps are searched.
    :return:
    """
    db.cursor.execute(
            """UPDATE snp
                SET %s = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (close_eQT, str(row['SNPChr']), int(row['SNPPos'])-region, 
            int(row['SNPPos'])+region))


def loop_eQTL(db, eQTL_df, ID_eQT, eQT, close_eQT, region):
    """
    Loop over the eQT dataframe
    :param db:  The database object
    :param eQTL_df:  The eQT file as data frame
    :param ID_eQT:  String. This is the name of the parameter from the database that is set to an ID 
                    when that snp exactly matches an eQT.
    :param eQT:  String. This is the name of the parameter in the database that will 
                 be set to true when a snp matches an eQT exactly
    :param close_eQT: String. This is the name of the parameter in the database that is set 
                      to true when a snp falls within a certain region (before or after) a known eQT
    :param region: Region (front or back) an eQT, within which snps are searched.
    :return:
    """
    for index, row in eQTL_df.iterrows():
        print(index)
        set_value(db, ID_eQT, eQT, row)
        set_value_close_to(db, row, close_eQT, region)
    # Add to database
    db.mydb_connection.commit()


def main():
    # Database file
    db_path='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/db_laatste_copy.db' #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"#"D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    # File with eQTLs
    eQTL_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/genes_eQTL_etc/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt'#"D:/Hanze_Groningen/STAGE/eQTL/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # eQTM_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    # Database class
    db = Database(db_path) #sys.argv[1]
    # Read file
    eQTL_df = pd.read_csv(eQTL_path, sep='\t')
    # eQTM_df = pd.read_csv(eQTL_path, sep='\t')
    region = 100
    # add_value(db, 'ID_eQTL', 'eQTL', 'close_eQTL')
    # Add to database
    db.mydb_connection.commit()
    # Call fill_eQTL
    loop_eQTL(db, eQTL_df, 'ID_eQTL', 'eQTL', 'close_eQTL', region)
    # Add to database
    db.mydb_connection.commit()
    # Close connection cursor and mydb_connection
    db.close()


if __name__ == '__main__':
    main()

