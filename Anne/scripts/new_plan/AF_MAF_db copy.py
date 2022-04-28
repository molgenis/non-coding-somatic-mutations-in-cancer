from Database import Database
import sys
import multiprocessing as mp
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config




def get_snps(db):
    print('GET SNPS')
    db.cursor.execute("""
                    SELECT COUNT(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    print('ID')
    for res in results:
        print(f'{res[0]} - {res}')

    db.cursor.execute("""
                    SELECT COUNT(DISTINCT ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    print('UNI ID')
    for res in results:
        print(f'{res[0]} - {res}')

    db.cursor.execute("""
                    SELECT count(*) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    print('COL')
    for res in results:
        print(f'{res[0]} - {res}')
        max_id = res[0]


    db.cursor.execute("""
                    SELECT MAX(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    print('MAX')
    for res in results:
            print(f'{res[0]} - {res}')

    db.cursor.execute("""
                    SELECT MIN(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    print('MIN')
    for res in results:
        print(f'{res[0]} - {res}')

    

    return max_id


def cal_AF(db):
    # COUNT donors
    db.cursor.execute("""
                    SELECT COUNT(DISTINCT ID)
                    FROM donor
                """ )
    count_donor = db.cursor.fetchall()
    for don in count_donor:
        c_donors = don[0]
    print('DONOR', c_donors)
    # max_id = get_snps(db)
    
    #https://stackoverflow.com/questions/30649873/how-do-i-count-distinct-combinations-of-column-values

    




def main():
    config = get_config()
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    cal_AF(db)
      



if __name__ == '__main__':
    main()