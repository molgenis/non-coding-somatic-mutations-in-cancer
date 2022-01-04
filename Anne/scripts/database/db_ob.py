import sqlite3
import glob
import pandas as pd
import sys
import io
import os


class Project(object):
    def __init__(self, ID, project_ID):
        self.ID = ID
        self.project_ID = project_ID

    def print_values(self):
        return f'{self.ID}\t{self.project_ID}'

class Donor(object):
    def __init__(self, ID, donor_ID, project_ID):
        self.ID = ID
        self.donor_ID = donor_ID
        self.project_ID = project_ID

    def print_values(self):
        return f'{self.ID}\t{self.donor_ID}\t{self.project_ID}'

class SNP(object):
    def __init__(self, ID, chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                in_transcript, in_coding, in_exon, platform, seq_strategy, tissue_id):
        self.ID = ID
        self.chr = chr
        self.pos_start = pos_start
        self.pos_end = pos_end
        self.ref = ref
        self.alt = alt
        self.genome_version = genome_version
        self.depth = depth
        self.in_transcript = in_transcript
        self.in_coding = in_coding
        self.in_exon = in_exon
        self.platform = platform
        self.seq_strategy = seq_strategy
        self.tissue_id = tissue_id

    def print_values(self):
        return f'{self.ID}\t{self.chr}\t{self.pos_start}\t{self.ref}\t{self.alt}'


class DonorSNP(object):
    def __init__(self, donor_ID, donor_project_ID, snp_ID):
        self.donor_ID = donor_ID
        self.donor_project_ID = donor_project_ID
        self.snp_ID = snp_ID

    def print_values(self):
        return f'{self.donor_ID}\t{self.donor_project_ID}\t{self.snp_ID}'

class Database:
    """
    
    """

    def __init__(self, db_name='D:/Hanze_Groningen/STAGE/TEST_DEL/db_test'):
        """
        Constructor
        :param db_name: 
        """
        self.name = db_name
        self.mydb_connection = self.connect()
        self.mydb_connection.row_factory = sqlite3.Row
        self.cursor = self.mydb_connection.cursor()

    def connect(self):
        """
        
        :return: 
        """
        try:
            return sqlite3.connect(self.name)
        except sqlite3.Error as er:
            print("Error while connecting to sqlite", er)
            pass

    def close(self):  # __del__
        """
        
        :return: 
        """
        self.cursor.close()
        self.mydb_connection.close()

    def check_table(self, table, column='ID'):
        """
        
        :return: 
        """
        print('value check')
        for row in self.cursor.execute(f"""SELECT * 
                                        FROM {table}
                                        ORDER BY {column} ASC"""):
            if table == 'project':
                print(Project(*row).print_values())
            elif table == 'donor':
                print(Donor(*row).print_values())
            elif table == 'snp':
                print(SNP(*row).print_values())
            elif table == 'donor_has_snp':
                print(DonorSNP(*row).print_values())

    def count_values(self, column, table):
        """

        :param column:
        :param table:
        :return:
        """
        self.cursor.execute(f"""
                            SELECT {column}, COUNT(*)
                            FROM {table}
                            GROUP BY {column};
                            """)
        results = self.cursor.fetchall()
        print(f'---{column}')
        for res in results:
            print(f'{res[0]} - {res[1]}')


    def create_db(self):
        """

        :return:
        """
        # Create tables
        self.cursor.execute("""
            -- -----------------------------------------------------
            -- Table `project`
            -- -----------------------------------------------------
            CREATE TABLE IF NOT EXISTS `project`(
                `ID` INTEGER PRIMARY KEY AUTOINCREMENT,
                `project_ID` VARCHAR(45) NOT NULL  
            )
            """)
        self.cursor.execute("""
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
        self.cursor.execute("""
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
        self.cursor.execute("""
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
        # return self.cursor

    def fill_database(self, df):
        """

        :param df:
        :return:
        """
        # Loop over set of project_ids and add it to the database
        for project_id in list(set(df['project_id'])):
            print(project_id)
            self.cursor.execute("""INSERT INTO project (project_ID) 
                            VALUES ('%s')""" % (str(project_id)))
            # Committing the current transactions
            self.mydb_connection.commit()
            # Get the last ID (private key of the project table) used
            last_id_project = self.cursor.lastrowid
            print(f"{last_id_project} - {project_id}")
            # Filter dataframe on project_id
            select_project = df.loc[df['project_id'] == project_id]
            print(f"donors: {len(set(select_project['donor_id']))}")
            # Loop over set of donor_ids in (last) project_id and add it to the database
            for donor_id in list(set(select_project['donor_id'])):
                self.cursor.execute("""INSERT INTO donor (donor_ID, project_ID)
                                VALUES ('%s', %s)""" % (str(donor_id), int(last_id_project)))
                # Committing the current transactions
                self.mydb_connection.commit()
                # Get the last ID (private key of the donor table) used
                last_id_donor = self.cursor.lastrowid
                # Filter dataframe on donor_id
                select_donor = select_project[select_project['donor_id'] == donor_id]
                # Loop over rows in dataframe (select_donor)
                for index, row in select_donor.iterrows():
                    # See if an SNP already exists with these values
                    self.cursor.execute(
                        """SELECT *
                        FROM snp
                        WHERE chr = '%s' AND pos_start = %s AND pos_end = %s AND 
                        ref = '%s' AND alt = '%s' AND genome_version = '%s' 
                        AND platform = '%s' AND seq_strategy = '%s';""" %
                        (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                         str(row['ref']), str(row['alt']), str(row['genome_version']),
                         str(row['platform']), str(row['seq_strategy'])))
                    check_snp = self.cursor.fetchall()
                    # If the SNP does not exist add it to the database
                    if not check_snp:
                        self.cursor.execute("""
                            INSERT INTO snp (chr, pos_start, pos_end, ref, alt, genome_version, depth, 
                                        platform, seq_strategy, tissue_id)
                            VALUES ('%s', %s, %s, '%s', '%s', '%s', %s, '%s', '%s', '%s')""" %
                                            (str(row['chr']), int(row['pos_start']), int(row['pos_end']),
                                             str(row['ref']), str(row['alt']), str(row['genome_version']),
                                             int(row['depth']),
                                             str(row['platform']), str(row['seq_strategy']), str(row['tissue_id'])))
                        # Get the last ID (private ket of the snp table) used
                        last_id_snp = self.cursor.lastrowid
                        # Fill the table donor_has_snp
                        self.cursor.execute("""
                            INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                            VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(last_id_snp)))
                        # Committing the current transactions
                        self.mydb_connection.commit()
                    # If the snp already exists insert the link between the donor and the snp by filling in
                    # the donor_has_snp table
                    else:
                        # Loop over the snp(s) it corresponds to, if all is well this is always 1 snp.
                        for info in check_snp:
                            # Get ID of the snp
                            id_snp = int(info['ID'])
                            # Check whether the combination donor and snp already exists
                            self.cursor.execute(
                                """SELECT *
                                FROM donor_has_snp
                                WHERE donor_project_ID = %s AND donor_ID = %s AND snp_ID = %s;""" %
                                (int(last_id_project), int(last_id_donor), int(id_snp))
                            )
                            check_donor_snp = self.cursor.fetchall()
                            # If the combination donor and snp does not yet exist, fill in the table donor_has_snp
                            # with this combination
                            if not check_donor_snp:
                                self.cursor.execute("""
                                    INSERT INTO donor_has_snp (donor_project_ID, donor_ID, snp_ID)
                                    VALUES (%s, %s, %s)""" % (int(last_id_project), int(last_id_donor), int(id_snp)))
                            # Committing the current transactions
                            self.mydb_connection.commit()

    def check_gene(self, path_fgene):
        """

        :param path_fgene:
        :return:
        """
        print('check_gene')
        gene_df = pd.read_csv(path_fgene, sep='\t')
        for index, row in gene_df.iterrows():
            print(index)
            self.cursor.execute(
                """UPDATE snp
                    SET in_transcript = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
            self.mydb_connection.commit()
            self.cursor.execute(
                """UPDATE snp
                    SET in_coding = TRUE
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
            self.mydb_connection.commit()
            exon_start = row['exonStarts'].rstrip(',').split(',')
            exon_end = row['exonEnds'].rstrip(',').split(',')
            print(f"COUNT - {row['exonCount']}")
            for i in range(int(row['exonCount'])):
                self.cursor.execute(
                    """UPDATE snp
                        SET in_exon = TRUE
                        WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
                    (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
                self.mydb_connection.commit()

    def compare_db(self):
        print()


def read_file(path, db):
    """
    
    :param path: 
    :param db: 
    :return: 
    """
    df = pd.read_csv(path, sep='\t')
    # Drop all duplicates (only depth and tissue_id may differ from all columns)
    df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
    db.fill_database(df)

def read_externDB(path, db):
    dbsnp = pd.read_csv(path, sep='\t')
    print(dbsnp.head())
    # db.compare_db()



def main():
    """
    
    :return: 
    """
    path = "D:/Hanze_Groningen/STAGE/db/files/ALL-US_db.tsv"
    path_fgene = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    dbSNP_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    db = Database()
    # db.create_db()
    # read_file(path, db)
    # db.check_gene(path_fgene)
    # db.check_donor()
    # db.count_values('in_transcript', 'snp')
    # db.count_values('in_coding', 'snp')
    # db.count_values('in_exon', 'snp')
    # db.check_table('snp')
    read_externDB(dbSNP_path, db)
    

    db.close()


main()