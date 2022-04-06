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



def cal_AF(db, type_GT):
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

    for ID in range(1, 10000): #TODO max id
        print(ID)
        # COUNT GT (but not if none/null)
        db.cursor.execute("""
                        SELECT %s, snp_ID, COUNT(*) as used_count
                        FROM donor_has_snp
                        WHERE snp_ID = %s AND (%s = 0 OR %s = 1 OR %s = 2)
                        GROUP BY %s, snp_ID
                        ORDER BY snp_ID;
                    """ % (type_GT, int(ID), type_GT, type_GT, type_GT, type_GT))
        results = db.cursor.fetchall()
        dict_count = dict()
        donor_count = 0
        donor_count_uniek = 0
        for res in results:
            if res[type_GT] in dict_count:
                dict_count[res[type_GT]] = dict_count[res[type_GT]] + 1
            else:
                dict_count[res[type_GT]] = 1
           
            # COUNT donors
            db.cursor.execute("""
                            SELECT COUNT(donor_ID)
                            FROM donor_has_snp 
                            WHERE snp_ID = %s AND (%s = 0 OR %s = 1 OR %s = 2);
                        """ %
            (int(res['snp_ID']), type_GT, type_GT, type_GT))
            count_donor = db.cursor.fetchall()
            for don in count_donor:
                donor_count += don[0]
                # print(f'{don[0]} - {don}')

            db.cursor.execute("""
                            SELECT COUNT(DISTINCT donor_ID)
                            FROM donor_has_snp 
                            WHERE snp_ID = %s AND (%s = 0 OR %s = 1 OR %s = 2);
                        """ %
            (int(res['snp_ID']), type_GT, type_GT, type_GT))
            count_donor = db.cursor.fetchall()
            for don in count_donor:
                donor_count_uniek += don[0]
                # print(f'{don[0]} - {don}')
        print('dict count', dict_count, 'donor count uniek', donor_count_uniek)
        ALT_all = 0
        for key, value in dict_count.items():
            ALT_all += (int(key) * int(value))
        AF = ALT_all / (donor_count_uniek * 2)
        print(f"AF: {AF} = {ALT_all} / ({donor_count_uniek} * 2)")

        AF_whole = ALT_all / (c_donors * 2)
        print(f"AF whole: {AF_whole} = {ALT_all} / ({c_donors} * 2)")
        # if donor_count != donor_count_uniek:
        #     print(ID)
        #     print('donor count', donor_count)
        #     print('donor count uniek', donor_count_uniek)
        #     db.cursor.execute("""
        #                     SELECT *
        #                     FROM donor_has_snp
        #                     WHERE snp_ID = %s AND (GT = 0 OR GT = 1 OR GT = 2);
        #                 """ %
        #     (int(ID)))
        #     results = db.cursor.fetchall()
        #     for res in results:
        #         print(f"{res['GT']} - {res['snp_ID']} - {res['donor_ID']}")
        #     print('++++++++++++++++++++++++')
        

    

    
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
    cal_AF(db, 'GT')
      



if __name__ == '__main__':
    main()