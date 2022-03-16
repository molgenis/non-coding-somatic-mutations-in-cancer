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