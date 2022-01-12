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

    def count_values(self, column, table, where):
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
        print(f'FILTER: {where}')
        self.cursor.execute(f"""
                            SELECT {column}, COUNT(*)
                            FROM {table}
                            WHERE {where}
                            GROUP BY {column};
                            """)
        results = self.cursor.fetchall()
        for res in results:
            print(f'{res[0]} - {res[1]}')

        print('\n\n')

