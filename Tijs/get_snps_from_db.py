#!/usr/bin/env python3
"""
Query mySQL database of ICGC data for SNPs.

MIT License

Copyright (c) 2022 Tijs van Lieshout

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Uses:
<The terminal interactions with this script go here>
"""

# Metadata
__title__ = "Query mySQL database of ICGC data for SNPs" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-04"
__updated__ = "2022-05-04"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 0.1
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse

import pandas as pd
import time

from utilities import get_config
from Database import Database


def main(args):
  config = get_config()
  db = Database(config['db_path'])

  if args.tables:
    if "," in args.tables:
      for table in args.tables.split(","):
        show_column_names(db, table)
        print(get_header_of_table(db, table))
    else:
      show_column_names(db, args.tables)
      print(get_header_of_table(db, args.tables))

  all_snps_df = get_all_snps_from_db(db)
  print(all_snps_df.head())
  print("Exporting to temporary tsv-file")
  all_snps_df.to_csv("tmp.tsv", sep='\t')
  db.close()


def get_tables(db):
  table = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table';", db.mydb_connection)
  return table


def show_column_names(db, table_name):
  for result in get_column_names(db, table_name):
    print(result[:])


def get_column_names(db, table_name):
  db.cursor.execute(f"PRAGMA table_info({table_name});")
  return db.cursor.fetchall()


def get_header_of_table(db, table):
  table = pd.read_sql(f"SELECT * FROM {table} LIMIT 10;", db.mydb_connection)
  return table

def get_all_snps_from_db(db):
  """
  Original author: Anne van Ewijk
  """
  print("Retrieving SNPs from table... (this might take around 1 min.)")
  start = time.time()
  table = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                               sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, 
                               snp.chr, snp.pos_start, snp.pos_end 
                        FROM project, sum_dosage_GT, snp 
                        WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                              sum_dosage_GT.donor_project_ID = project.ID AND 
                              (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                              sum_dosage_GT.total_read_count_sum >= 33;''', db.mydb_connection)
  print(f"Retrieved SNPs from table in {time.time()-start:.0f} seconds!")
  return table

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-t", "--tables", type=str, required=False, help="Name of tables to query. Add multiple tables by , (e.g. snp,donor")
  #parser.add_argument("-b", "--example_bool", action="store_true", required=False, help="Help goes here") 
  args = parser.parse_args()
  
  main(args)
