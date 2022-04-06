#!/usr/bin/env python3
from os import sep
import pandas as pd
import sys
import math # nan
import numpy as np

from Database import Database


def main():
    # Make Database object
    db = Database('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydb_L.db') #00Database_internship_UPDATE2.0.db.db')
    projectset = set()
    db.cursor.execute("""SELECT *
                            FROM project""")
    print('yo1')
    projects = db.cursor.fetchall()
    for proj in projects:
        projectset.add(proj['project_ID'])
        db.cursor.execute("""
                    SELECT COUNT(DISTINCT donor_ID)
                    FROM 'donor_has_snp'
                    WHERE donor_project_ID = %s;
                    """ % (proj['ID']))
        results = db.cursor.fetchall()
        print('UNI ID')
        for res in results:
            print(f"{proj['project_ID']} ({proj['ID']}) - {res[0]} - {res}")
        
    print(f'len proj: {len(projectset)}')
    print('yo2')
    donorset = set()
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    for don in donors:
        donorset.add(don['donor_ID'])
    print(f'len donor: {len(donorset)}')

if __name__ == '__main__':
    main()
