from Database import Database

import pandas as pd

def check_number_donors_cancer(db):
    db.cursor.execute(
        """SELECT cancer
            FROM project;""" )    
    results = db.cursor.fetchall()
    cancer_set = set()
    for res in results:
        cancer_set.add(res['cancer'])

    for can in list(cancer_set):
        db.cursor.execute(
            """SELECT COUNT(DISTINCT donor.ID) 
                FROM project, donor
                WHERE project.cancer = '%s' and project.ID=donor.project_ID;""" %
            (str(can)))
        
        results = db.cursor.fetchall()
        for res in results:
            print(f'{can} - {res[0]}')


def check_number_donors_project(db):
    db.cursor.execute(
        """SELECT ID
            FROM project;""" )    
    results = db.cursor.fetchall()
    cancer_set = set()
    for res in results:
        cancer_set.add(res['ID'])

    for can in list(cancer_set):
        db.cursor.execute(
            """SELECT project.ID, project.project_ID, project.cancer, COUNT(DISTINCT donor.ID) 
                FROM project, donor
                WHERE project.ID = '%s' and project.ID=donor.project_ID;""" %
            (str(can)))
        
        results = db.cursor.fetchall()
        for res in results:
            print(f"{res[0]} - {res[1].replace('-', '_')} - {res[2]} - {res[3]}")


def main():
    #
    path_db = 'D:/Hanze_Groningen/STAGE/db_laatste_copy.db' #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    check_number_donors_cancer(db)
    check_number_donors_project(db)



    
    
    




if __name__ == '__main__':
    main()