from Database import Database
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from multiprocessing import Pool, Queue
import multiprocessing as mp
import math

import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



# def search_close_snp(db, df_strong_eQTL, region, config):
#     """
    
#     """
#     f = open(config['distance_eqtl_path'], 'w') #D:/Hanze_Groningen/STAGE/eQTL/
#     f.write(f"SNP_eQTL\chr\pos\Gene\tsnp_ID\tchr_snp\tpos_snp\tdistance\n")
#     for index, row in df_strong_eQTL.iterrows():
#         print('----------')
#         print(str(row['SNPChr']), int(row['SNPPos']), int(row['SNPPos']))
#         db.cursor.execute(
#         """SELECT ID, chr, pos_start, abs(CAST(%s AS REAL) - CAST(pos_start AS REAL)) AS distance_eQTL
#             FROM snp
#             WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
#             ORDER BY distance_eQTL
#             LIMIT 1;""" %
#         (int(row['SNPPos']), str(row['SNPChr']), int(row['SNPPos']-region), int(row['SNPPos']+region)))
    
#         results = db.cursor.fetchall()

#         # for res in results:
#         #     print(f"{res[0]}-{res[1]}-{res[2]}-{res[3]}")

#         if len(results) != 0:
#             for res in results:
#                 f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\n")
#         else:
#             print('MILJOEN')
#             db.cursor.execute(
#             """SELECT ID, chr, pos_start, abs(CAST(%s AS REAL) - CAST(pos_start AS REAL)) AS distance_eQTL
#                 FROM snp
#                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
#                 ORDER BY distance_eQTL
#                 LIMIT 1;""" %
#             (int(row['SNPPos']), str(row['SNPChr']), int(row['SNPPos']-1000000), int(row['SNPPos']+1000000)))
        
#             results = db.cursor.fetchall()
#             if len(results) != 0:
#                 for res in results:
#                     f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\n")
#             else:
#                 f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t-\t-\t-\t-\n")
#     f.close()

def multiprocess_search_close_snp(df_strong_eQTL, region, path_save_file, config):
    """
    
    """
    path_db = config['database_2'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" 
    # Database connection
    db = Database(path_db)
    f = open(path_save_file, 'w') #D:/Hanze_Groningen/STAGE/eQTL/
    f.write(f"SNP_eQTL\tchr\tpos\tGene\tsnp_ID\tchr_snp\tpos_snp\tdistance\n")
    for index, row in df_strong_eQTL.iterrows():
        print('----------')
        print(str(row['SNPChr']), int(row['SNPPos']), int(row['SNPPos']))
        db.cursor.execute(
        """SELECT ID, chr, pos_start, abs(CAST(%s AS REAL) - CAST(pos_start AS REAL)) AS distance_eQTL
            FROM snp
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
            ORDER BY distance_eQTL
            LIMIT 1;""" %
        (int(row['SNPPos']), str(row['SNPChr']), int(row['SNPPos']-region), int(row['SNPPos']+region)))
    
        results = db.cursor.fetchall()

        # for res in results:
        #     print(f"{res[0]}-{res[1]}-{res[2]}-{res[3]}")

        if len(results) != 0:
            for res in results:
                f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\n")
        else:
            print('MILJOEN')
            db.cursor.execute(
            """SELECT ID, chr, pos_start, abs(CAST(%s AS REAL) - CAST(pos_start AS REAL)) AS distance_eQTL
                FROM snp
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                ORDER BY distance_eQTL
                LIMIT 1;""" %
            (int(row['SNPPos']), str(row['SNPChr']), int(row['SNPPos']-1000000), int(row['SNPPos']+1000000)))
        
            results = db.cursor.fetchall()
            if len(results) != 0:
                for res in results:
                    f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\n")
            else:
                f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t-\t-\t-\t-\n")
    f.close()

def main():
    config = get_config()
    # #
    # path_db = config['database'] #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # # Database connection
    # db = Database(path_db)
    path_strong_eQTL = config['strong_eqtl_path']
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t')
    region = 10000

    cpus = mp.cpu_count()

    list_paths_save = list()
    arg_multi_list = []
    for i in range(0, len(df_strong_eQTL)+1000, 1000):
        df = df_strong_eQTL.iloc[i:i+1000]
        path_save_file = f'{config["genes_eQTL_etc"]}distance_eqtl_{i}.tsv'
        list_paths_save.append(path_save_file)
        arg_multi_list.append((df, region, path_save_file, config))

    pool = Pool(processes=cpus)
    pool.starmap(func=multiprocess_search_close_snp, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # search_close_snp(db, df_strong_eQTL, region, config)


if __name__ == '__main__':
    main()