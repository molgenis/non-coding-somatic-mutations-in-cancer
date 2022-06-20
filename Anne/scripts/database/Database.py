#!/usr/bin/env python3
import sqlite3
# import matplotlib.pyplot as plt
import numpy as np


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
        print(self.mydb_connection)
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

    def count_values(self, column, table, where, title):
        """
        Counts rows in the column that possibly meet a given 'where' statement in the database
        :param column: The selected column of the table in the database
        :param table:  The selected table in the database
        :param where:  The 'where' statement in the SQL code
        :return:
        """
        # Count the rows without the 'where' statement
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
        # Count the rows with the 'where' statement
        self.cursor.execute(f"""
                            SELECT {column}, COUNT(*)
                            FROM {table}
                            WHERE {where}
                            GROUP BY {column};
                            """)
        results = self.cursor.fetchall()
        for res in results:
            print(f'{res[0]} - {res[1]}')
        # self.make_plot_true_false(results, title)
        return results
        

    def make_plot_true_false(self, results, title):
        if results[0][0] == 0:
            names = ['False', 'True']
        elif results[0][0] == 1:
            names = ['True', 'False']
        
        if len(results) == 2:
            values = [results[0][1], results[1][1]]
        elif len(results) == 1:
            values = [results[0][1], 0]


        # plt.bar(names, values)
        # # Add title and axis names
        # plt.title(title)
        # plt.xlabel('Categories')
        # plt.ylabel('Counts')
        #
        # plt.show()

        

        # print(results[1][1])
        # for index, res in enumerate(results):
        #     if res[0] == 1:
        #         #TRUE
        #         plt.bar(res[0], res[1], color = 'g')
        #     else:
        #         plt.bar(res[0], res[1], color = 'r')
        
        # plt.show()
