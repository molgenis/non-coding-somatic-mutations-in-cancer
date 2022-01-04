import sqlite3
import glob
import pandas as pd
import sys

try:
    # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
    mydb_connection = sqlite3.connect("D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene.db")
    # Setting row_factory property of connection object to
    # sqlite3.Row(sqlite3.Row is an implementation of row_factory)
    mydb_connection.row_factory = sqlite3.Row
    # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
    cursor = mydb_connection.cursor()
    

    cursor.execute('SELECT * FROM donor')
    for x in cursor:
        print(f"{x['donor_ID']} - {x['project_ID']} ")
except sqlite3.Error as er:
    print("Error while connecting to sqlite", er)
finally:
    if mydb_connection:
        mydb_connection.close()
        print("The SQLite connection is closed")