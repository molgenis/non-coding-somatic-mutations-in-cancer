#!/usr/bin/env python3
from os import sep
import pandas as pd
import sys
import math # nan
import numpy as np

from Database import Database


def main():
    # Make Database object
    db = Database('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/00Database_internship_UPDATE2.0.db.db')
    projectset = set()
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    for proj in projects:
        print(proj['project_ID'])
        projectset.add(proj['project_ID'])
    print(f'len proj: {len(projectset)}')

    donorset = set()
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    for don in donors:
        donorset.add(don['donor_ID'])
    print(f'len proj: {len(donorset)}')

if __name__ == '__main__':
    main()
