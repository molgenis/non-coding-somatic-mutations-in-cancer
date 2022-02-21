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
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    for proj in projects:
        print(proj['project_ID'])

if __name__ == '__main__':
    main()
