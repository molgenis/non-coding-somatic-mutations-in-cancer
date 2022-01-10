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
        print(self.mydb_connection)
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




# def read_file(path, db):
#     """
    
#     :param path: 
#     :param db: 
#     :return: 
#     """
#     df = pd.read_csv(path, sep='\t')
#     # Drop all duplicates (only depth and tissue_id may differ from all columns)
#     df = df.drop_duplicates(subset=df.columns.difference(['depth', 'tissue_id']))
#     db.fill_database(df)

# def read_externDB(path, db):
#     dbsnp = pd.read_csv(path, sep='\t')
#     print(dbsnp.head())
#     # db.compare_db()



# def main():
#     """
    
#     :return: 
#     """
#     path = "D:/Hanze_Groningen/STAGE/db/files/ALL-US_db.tsv"
#     path_fgene = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
#     dbSNP_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
#     db = Database()
#     # db.create_db()
#     # read_file(path, db)
#     # db.check_gene(path_fgene)
#     # db.check_donor()
#     # db.count_values('in_transcript', 'snp')
#     # db.count_values('in_coding', 'snp')
#     # db.count_values('in_exon', 'snp')
#     # db.check_table('snp')
#     # read_externDB(dbSNP_path, db)
    

#     db.close()


# main()