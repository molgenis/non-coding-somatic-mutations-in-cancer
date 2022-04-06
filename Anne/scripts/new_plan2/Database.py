#!/usr/bin/env python3
import sqlite3


class Database:
    """
    
    """

    def __init__(self, db_name):
        """
        Constructor
        :param db_name: The name of the database
        """
        self.name = db_name
        self.mydb_connection = self.connect()
        # print(self.mydb_connection)
        self.mydb_connection.row_factory = sqlite3.Row
        self.cursor = self.mydb_connection.cursor()

    def connect(self):
        """
        Connect to the database
        :return: 
        """
        try:
            return sqlite3.connect(self.name)
        except sqlite3.Error as er:
            print("Error while connecting to sqlite", er)
            pass

    def close(self):  # __del__
        """
        Close the connections of the cursor and mydb_connection
        :return: 
        """
        self.cursor.close()
        self.mydb_connection.close()

    def get_projects(self):
        """
        Searches for all projects in the database and converts them into a specific dictionary. (key: ID, value: cancer)
        :param db:             The database object
        :return: project_dict: Dictionary with as key project ID (automatically generated) and as value the type of cancer
                               that belongs to it.
        """
        # Select all projects in the database
        self.cursor.execute("""SELECT *
                                FROM project""")
        projects = self.cursor.fetchall()
        project_dict = dict()
        for proj in projects:
            # Add to project_dict as key the auto-generated ID and as value the type of cancer that belongs to that project
            project_dict[proj['ID']] = proj['cancer']
        return project_dict


    def get_donors(self, project_dict):
        """
        Search for all donors in the database and make different things with them (lists and dictionaries).
        :param db:                  The database object
        :param project_dict:        Dictionary with as key project ID (automatically generated) and as value the type of
                                    cancer that belongs to it.
        :return: donor_list:        List of donor names (to be used later as rows in the sparse matrix)
                donor_dict:         A dictionary with as key the automatically generated donor ID and as value the donor
                                    IDs that are used in the research.
                donor_cancer_list:  List of cancers. This list has the same order as donor_list.
        """
        # Select all donors in the database
        self.cursor.execute("""SELECT *
                                FROM donor""")
        donors = self.cursor.fetchall()
        donor_dict = dict()
        donor_list = list()
        donor_cancer_list = list()
        for donor in donors:
            # Add to donor_dict as key the auto-generated ID and as value the donor ID as used in the study
            donor_dict[donor['ID']] = donor['donor_ID']
            donor_list.append(donor['donor_ID'])
            donor_cancer_list.append(project_dict[donor['project_ID']])
        return donor_list, donor_dict, donor_cancer_list

    def get_snps(self, chr, start_pos, end_pos):
        """
        Finds the snps between certain positions and puts the snp IDs in a list.
        :param chr:         Chromosoom number or letter
        :param start_pos:   The start position
        :param end_pos:     The end position
        :return:
        """
        # Find all snps that are on a certain chromosome in a certain region (between a certain start and stop)
        self.cursor.execute("""
                        SELECT ID
                        FROM 'snp'
                        WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                        """ %
                        (str(chr), int(start_pos), int(end_pos)))
        results_before = self.cursor.fetchall()
        # Make snp_id_list
        snp_id_list = list()
        for res_bef in results_before:
            # Add ID to snp_id_list
            snp_id_list.append(res_bef['ID'])
        return snp_id_list

