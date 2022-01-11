import sqlite3
import glob
import pandas as pd
import sys
import io
import os

from db_ob import Database

db = Database(sys.argv[1])

mydb_connection = db.mydb_connection
cursor = db.cursor

f = open(f"{sys.argv[2]}chr{sys.argv[3]}_db.vcf", "a")
f.write("##fileformat=VCFv4.0\n")
f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")

cursor.execute(
            """SELECT *
            FROM snp
            WHERE seq_strategy = 'WGS' AND chr = '%s' ;""" % #AND pos_start = %s AND pos_end = %s AND 
            # ref = '%s' AND alt = '%s' 
            (str(sys.argv[3]))) #, int(sys.argv[2]), int(sys.argv[3]),
            #str(sys.argv[4]), str(sys.argv[5])))
results = cursor.fetchall()
for res in results:
    f.write(f"chr{res['chr']}\t{res['pos_start']}\t.\t{res['ref']}\t{res['alt']}\t.\t.\t.\t.\n")

f.close()
db.close()
