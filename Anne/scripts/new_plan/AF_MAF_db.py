from Database import Database
import sys
import multiprocessing as mp


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

    

    # # SNP
    # db.cursor.execute("""
    #                 SELECT ID
    #                 FROM 'snp'
    #                 """)
    # snps = db.cursor.fetchall()
    # snp_list = list()
    # for snp in snps:
    #     snp_list.append(snp['ID'])
    # print(len(snp_list))
    return max_id


# def multi_pro(db):
#     print('MULTI')
#     print(f'len: {len(snp_list)}')
#     print(f'len set: {len(set(snp_list))}')
#     cpus = mp.cpu_count()
#     # round up
#     n = math.ceil(len(snp_list) / cpus)
#     print(f'CPUS {cpus}')
#     print(n)
#     arg_multi_list = []
#     for i in range(0, len(snp_list), n):
#         name_file = f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/UMAP/chr{key}_{i}.tsv.gz' #f'D:/Hanze_Groningen/STAGE/NEW PLAN/chr{key}_{i}.tsv.gz' #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/UMAP/
#         arg_multi_list.append((snp_list[i: i+n], donor_dict, donor_list, step_index, sparseMatrix, step_list, donor_cancer_list, name_file))

#     pool = Pool(processes=cpus)
#     pool.starmap(func=get_number_snps_region, iterable=arg_multi_list)
#     pool.close()
#     pool.join()



def cal_AF(db):
    # max_id = get_snps(db)
    # for id in range(1, max_id):
    #     print()
    # print('hoi')

    # db.cursor.execute("""
    #                 SELECT GT, COUNT(*) as used_count
    #                 FROM donor_has_snp
    #                 GROUP BY GT;
    #             """) # %
    # # (int(2222))) #WHERE snp_ID = %s
    # print('hoi2')
    # results = db.cursor.fetchall()
    # print('start')
    # for res in results:
    #     print(f"{res[0]}-{res[1]}")
    #https://stackoverflow.com/questions/30649873/how-do-i-count-distinct-combinations-of-column-values

    db.cursor.execute("""
                    SELECT *
                    FROM donor_has_snp
                    WHERE snp_ID = %s AND GT = %s;
                """ %
    (int(259183), 1))
    results = db.cursor.fetchall()
    for res in results:
        print(f"{res['GT']} - {res['snp_ID']} - {res['donor_ID']}")

    print('######################')

    db.cursor.execute("""
                    SELECT GT, snp_ID, COUNT(*) as used_count
                    FROM donor_has_snp
                    WHERE snp_ID = %s
                    GROUP BY GT, snp_ID
                    ORDER BY snp_ID;
                """ %
    (int(259183)))
    results = db.cursor.fetchall()

    for res in results:
        print('---------------')
        print(f"{res['GT']} - {res['snp_ID']} - {res['used_count']}")
        if isinstance(res['GT'], type(None)):
            GT = 'NULL'
        else:
            GT = res['GT']
        # COUNT donors
        print('DONOR')
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID, COUNT(*) as donor_c
                        FROM donor_has_snp
                        WHERE snp_ID = %s AND GT = %s
                        GROUP BY snp_ID
                        ORDER BY snp_ID;
                    """ %
        (int(res['snp_ID']), GT))
        count_donor = db.cursor.fetchall()
        for cd in count_donor:
            print(f"{cd['snp_ID']} - {cd['donor_c']}")
        print('UNIEK')
        db.cursor.execute("""
                        SELECT COUNT(DISTINCT donor_ID)
                        FROM donor_has_snp 
                        WHERE snp_ID = %s AND GT = %s;
                    """ %
        (int(res['snp_ID']), GT))
        count_donor = db.cursor.fetchall()
        for cd in count_donor:
            print(f'{cd[0]} - {cd}')
        print('NIET UNIEK')
        db.cursor.execute("""
                        SELECT COUNT(donor_ID)
                        FROM donor_has_snp
                        WHERE snp_ID = %s AND GT = %s;
                    """ %
        (int(res['snp_ID']), GT))
        count_donor = db.cursor.fetchall()
        for cd in count_donor:
            print(f'{cd[0]} - {cd}')

    

    
    # db.cursor.execute("""
    #                 SELECT GT, snp_ID
    #                 FROM donor_has_snp
    #                 WHERE snp_ID = 14868
    #             """)
    # results = db.cursor.fetchall()
    # dict_GT = dict()
    # for res in results:
    #     print(res['GT'])
    #     if res['GT'] in dict_GT:
    #         dict_GT[res['GT']] = dict_GT[res['GT']] + 1
    #     else:
    #         dict_GT[res['GT']] = 1

    # print(dict_GT)
    

    # for index, ID_snp in enumerate(snp_list):
    #     if index % 100 == 0:
    #         print(index)
    #     # db.cursor.execute("""
    #     #                 SELECT GT
    #     #                 FROM donor_has_snp
    #     #                 WHERE snp_ID = '%s'
    #     #             """ %
    #     # (int(ID_snp)))
    #     # donor_had_snps = db.cursor.fetchall()
    #     # GT_dict = dict()
    #     # for GT in donor_had_snps:
    #     #     if GT['GT'] in GT_dict:
    #     #         GT['GT'] = GT['GT'] + 1
    #     #     else:
    #     #         GT['GT'] = 1
    #     db.cursor.execute("""
    #                     SELECT
    #                         GT,
    #                         COUNT(GT)
    #                     FROM
    #                         donor_has_snp
    #                     GROUP BY
    #                         GT
    #                     WHERE snp_ID = '%s';
    #                 """ %
    #     (int(ID_snp)))
    #     results = db.cursor.fetchall()
    #     for res in results:
    #         print(f'{res[0]} - {res[1]}')






def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    cal_AF(db)
      



if __name__ == '__main__':
    main()