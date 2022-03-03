from Database import Database
import time
import pandas as pd
from multiprocessing import Pool, Queue
import multiprocessing as mp
# Python program to create
# sparse matrix using csr_matrix()
  
# Import required package
import numpy as np
import scipy.sparse
from scipy.sparse import csr_matrix


# def steps_snp(db, key, i_start, i, donor_dict, donor_list,
#                 sparseMatrix, step_index):
#     # SNP
#     db.cursor.execute("""
#                     SELECT ID
#                     FROM 'snp'
#                     WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
#                     """ %
#         (str(key), int(i_start), int(i)))
#     snps = db.cursor.fetchall()
#     snp_list = list()
#     for snp in snps:
#         snp_list.append(snp['ID'])
    
#     # donor_has_snp
#     for ID_snp in snp_list:
#         db.cursor.execute("""
#                         SELECT donor_ID
#                         FROM donor_has_snp
#                         WHERE snp_ID = '%s'
#                     """ %
#         (int(ID_snp)))
#         donor_had_snps = db.cursor.fetchall()
#         for d_h_s in donor_had_snps:
#             donor_ID = donor_dict[d_h_s['donor_ID']]
#             donor_index = donor_list.index(donor_ID)
#             sparseMatrix[donor_index, step_index] += 1
            
#     return sparseMatrix, donor_index


def get_projects(db):
    # PROJECT
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    project_dict = dict()
    for proj in projects:
        project_dict[proj['ID']] = proj['cancer']
    return project_dict



def get_donors(db, project_dict):
    # DONOR
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    donor_dict = dict()
    donor_list = list()
    donor_cancer_list = list()
    for donor in donors:
        donor_dict[donor['ID']] = donor['donor_ID']
        donor_list.append(donor['donor_ID'])
        donor_cancer_list.append(project_dict[donor['project_ID']])
    return donor_list, donor_dict, donor_cancer_list



def get_snps(db, key, i_start, i):
    # SNP
    db.cursor.execute("""
                    SELECT ID
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    """ %
        (str(key), int(i_start), int(i)))
    snps = db.cursor.fetchall()
    snp_list = list()
    for snp in snps:
        snp_list.append(snp['ID'])
    return snp_list



def get_number_snps_region(db, snp_list, donor_dict, donor_list, step_index, sparseMatrix, count):
    print(len(snp_list))
    # donor_has_snp
    # TODO multiprocessen
    for ID_snp in snp_list:
        count += 1
        print(count)
        db.cursor.execute("""
                        SELECT donor_ID
                        FROM donor_has_snp
                        WHERE snp_ID = '%s'
                    """ %
        (int(ID_snp)))
        donor_had_snps = db.cursor.fetchall()
        for d_h_s in donor_had_snps:
            donor_ID = donor_dict[d_h_s['donor_ID']]
            donor_index = donor_list.index(donor_ID)
            sparseMatrix[donor_index, step_index] += 1
            
    return count, sparseMatrix, donor_index


def make_table(db, chr_length, steps, donor_list, donor_dict, donor_cancer_list):
    for key, value in chr_length.items():
        count = 0 #TODO REMOVE
        name_file = 'D:/Hanze_Groningen/STAGE/NEW PLAN/chr{key}.tsv.gz'
        i_start=0
        step_list = []
        for i in range(0,value,steps):
            step_list.append(f'{key}:{i_start}-{i}')
            if value == (i + (value % steps)):
                last_i = (i + (value % steps))
                step_list.append(f'{key}:{(i + 1)}-{last_i}')
            i_start = (i + 1)
        # Creating a len(donor_list) * len(step_list) sparse matrix
        sparseMatrix = csr_matrix((len(donor_list), len(step_list)), 
                                dtype = np.int8).toarray()
        i_start=0
        for i in range(0,value,steps):
            step_index = step_list.index(f'{key}:{i_start}-{i}')
            print(f'{i_start} - {i}')
            if i != 0:
                snp_list = get_snps(db, key, i_start, i)
                # TODO multprocessing
                # https://stackoverflow.com/questions/47900922/split-list-into-n-lists-and-assign-each-list-to-a-worker-in-multithreading
                q = Queue(snp_list)
                cpu = mp.cpu_count()-2
                for i in range(cpu):
                    pool = Pool(processes=cpu)
                    pool.starmap(func=make_table, iterable=(q,))

                # pool = Pool(processes=cpu)
                # pool.starmap(func=make_table, iterable=arg_multi_list)
                count, sparseMatrix, donor_index = get_number_snps_region(db, snp_list, donor_dict, donor_list, step_index, sparseMatrix, count)
                # sparseMatrix, donor_index = steps_snp(db, key, i_start, i, donor_dict, donor_list, sparseMatrix, step_index)
            if value == (i + (value % steps)):
                print('YOOOOOOOOOOOO')
                last_i = (i + (value % steps))
                print(f'---- {i_start} - {i + 1} - {last_i}')
                snp_list = get_snps(db, key, (i+1), last_i)
                count, sparseMatrix, donor_index = get_number_snps_region(db, snp_list, donor_dict, donor_list, step_index, sparseMatrix, count)
                # sparseMatrix, donor_index = steps_snp(db, key, (i + 1), last_i, donor_dict, donor_list, sparseMatrix, step_index)

            i_start = (i + 1)


        df = pd.DataFrame(data=sparseMatrix, columns=step_list)
        df['donor_id'] = donor_index
        df['cancer'] = donor_cancer_list
        df.set_index('donor_id')
        df.to_csv(name_file, sep="\t", index=True, encoding='utf-8', 
                    compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})




def main():
    start_time = time.time()
    # The steps for the region
    steps= 10000000
    #https://en.wikipedia.org/wiki/Human_genome
    chr_length = {'1':45000680}#{'1':248956422, '2':242193529} #, '3':198295559, '4':190214555,
                # '5':181538259, '6':170805979, '7':159345973, '8':145138636,
                # '9':138394717, '10':133797422, '11':135086622, '12':133275309,
                # '13':114364328, '14':107043718, '15':101991189, '16':90338345,
                # '17':83257441, '18':80373285, '19':58617616, '20':64444167,
                # '21':46709983, '22':50818468, 'X':156040895, 'Y':57227415} #{'17':83257441, '18':80373285, '19':58617616, '20':64444167}
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"
    # Database connection
    db = Database(path_db)
    project_dict = get_projects(db)
    donor_list, donor_dict, donor_cancer_list = get_donors(db, project_dict)    
    make_table(db, chr_length, steps, donor_list, donor_dict, donor_cancer_list)  

    



if __name__ == '__main__':
    main()