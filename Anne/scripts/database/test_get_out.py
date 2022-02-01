#!/usr/bin/env python3
import sqlite3
import pandas as pd
from Database import Database

# def main():
#     try:
#         # Returns a connection object that we will use to interact with the SQLite database held in the file test.db
#         mydb_connection = sqlite3.connect("D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene.db")
#         # Setting row_factory property of connection object to
#         # sqlite3.Row(sqlite3.Row is an implementation of row_factory)
#         mydb_connection.row_factory = sqlite3.Row
#         # Cursor object allow us to send SQL statements to a SQLite database using cursor.execute()
#         cursor = mydb_connection.cursor()
        

#         cursor.execute('SELECT * FROM project')
#         for x in cursor:
#             print(f"{x['ID']} - {x['project_ID']} ")
#     except sqlite3.Error as er:
#         print("Error while connecting to sqlite", er)
#     finally:
#         if mydb_connection:
#             mydb_connection.close()
#             print("The SQLite connection is closed")

# if __name__ == '__main__':
#     main()


db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
db = Database(db_path)
db.cursor.execute(f"""
                    SELECT in_coding, COUNT(*)
                    FROM 'snp'
                    WHERE seq_strategy = 'WGS'
                    """)
results = db.cursor.fetchall()
for res in results:
    print(f'{res[0]} - {res[1]}')