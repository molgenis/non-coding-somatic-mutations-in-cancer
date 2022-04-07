from Database import Database
import sys
import multiprocessing as mp


def get_snps(db):
    """
    Finds the highest and lowest snp ID.
    :param db:  The database object
    :return: max_snp_id: The highest snp id
             min_snp_id: The lowest snp id
    """
    print('GET SNPS')
    # MAX
    db.cursor.execute("""
                    SELECT MAX(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        max_snp_id = res[0]
    # MIN
    db.cursor.execute("""
                    SELECT MIN(ID) 
                    FROM snp;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        min_snp_id = res[0]
    return max_snp_id, min_snp_id


def add_value(db):
    """
    Adds values (AF and AF_whole) to the database (table snp).
    :param db:  The database object
    :return:
    """ #TODO explain AF and AF_whole
    # Add in_transcript
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `AF` BOOLEAN DEFAULT(FALSE)
                    """)
    # Add in_coding
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `AF_whole` BOOLEAN DEFAULT(FALSE)
                    """)
    # Committing the current transactions
    db.mydb_connection.commit()



def find_number_donors(db):
    """
    Finds the number of donors in the database
    :param db:  The database object
    :return:  c_donors: Number of donors in the database
    """
    # COUNT donors
    db.cursor.execute("""
                    SELECT COUNT(DISTINCT ID)
                    FROM donor
                """ )
    count_donor = db.cursor.fetchall()
    for don in count_donor:
        c_donors = don[0]
    print('DONOR', c_donors)
    return c_donors



def cal_AF(db, type_GT):
    """
    
    :param db:  The database object
    :param type_GT:
    :return:
    """
    c_donors = find_number_donors(db)
    max_snp_id, min_snp_id = get_snps(db)
    print(max_snp_id)
    
    #https://stackoverflow.com/questions/30649873/how-do-i-count-distinct-combinations-of-column-values

    for ID in range(1, 200): #TODO min_snp_id, max_snp_id
        print(ID)
        # COUNT GT (but not if none/null)
        db.cursor.execute("""
                        SELECT %s, snp_ID, donor_ID, COUNT(*) as used_count
                        FROM sum_dosage_GT
                        WHERE snp_ID = %s AND (%s = 0 OR %s = 1 OR %s = 2)
                        GROUP BY %s, snp_ID
                        ORDER BY snp_ID;
                    """ % (type_GT, int(ID), type_GT, type_GT, type_GT, type_GT))
        results = db.cursor.fetchall()
        snp_count_dict = dict()
        donor_count_dict = dict()
        donor_count = 0
        donor_count_uniek = 0
        for res in results:
            # GT
            if res[type_GT] in snp_count_dict:
                snp_count_dict[res[type_GT]] = snp_count_dict[res[type_GT]] + 1
            else:
                snp_count_dict[res[type_GT]] = 1
            # donor_ID
            if res['donor_ID'] in donor_count_dict:
                print('BESTAAT AL')
                donor_count_dict[res['donor_ID']] = donor_count_dict[res['donor_ID']] + 1
                donor_count += 1
            else:
                donor_count_dict[res['donor_ID']] = 1
                donor_count += 1
                donor_count_uniek += 1

        # if donor_count != donor_count_uniek:
        #     print()
            # print(snp_count_dict)
            # print('donor_count', donor_count, 'donor_count_uniek', donor_count_uniek)
        ALT_all = 0
        for key, value in snp_count_dict.items():
            ALT_all += (int(key) * int(value)) ##########################################
        if donor_count_uniek != 0 or donor_count != 0:
            # # donor uniek
            AF2 = ALT_all / (donor_count_uniek * 2)
            # donor
            AF = ALT_all / (donor_count * 2)
            # all donors
            AF_whole = ALT_all / (c_donors * 2)
        else:
            print(ALT_all, donor_count, donor_count_uniek)
            AF = 0
            AF2 = 0
            AF_whole = 0

        if AF != AF2:
            print(snp_count_dict)
            print('donor_count', donor_count) #, 'donor_count_uniek', donor_count_uniek)
            # # donor uniek
            # print(f"AF uniek: {AF2} = {ALT_all} / ({donor_count_uniek} * 2)")
            # donor                
            print(f"AF: {AF} = {ALT_all} / ({donor_count} * 2)")
            # all donors            
            print(f"AF whole: {AF_whole} = {ALT_all} / ({c_donors} * 2)")

        db.cursor.execute(
            """UPDATE snp
                SET AF = %s, AF_whole = %s
                WHERE ID = %s;""" % (AF, AF_whole, ID))
    # Add to database
    db.mydb_connection.commit()


def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    #'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydb_L.db' #
    # Database connection
    db = Database(path_db)
    add_value(db)
    type_GT = 'GT'
    cal_AF(db, type_GT)
    print('#########################################')
    type_GT = 'GT2'
    cal_AF(db, type_GT)



    # #CHECK
    # db.cursor.execute("""
    #                 SELECT snp_ID, donor_ID, total_read_count, mutant_allele_read_count
    #                 FROM donor_has_snp
    #                 WHERE donor_ID = 614;
    #             """)
    # results = db.cursor.fetchall()
    # snp_count_dict = dict()
    # for res in results:
    #     # GT
    #     if res['snp_ID'] in snp_count_dict:
    #         snp_count_dict[res['snp_ID']] = snp_count_dict[res['snp_ID']] + 1
    #     else:
    #         snp_count_dict[res['snp_ID']] = 1
    # print(snp_count_dict)
    # #1436355
    # db.cursor.execute(f"""
    #                 SELECT snp_ID, donor_ID, total_read_count, mutant_allele_read_count, {type_GT}
    #                 FROM donor_has_snp
    #                 WHERE donor_ID = 614 AND snp_ID = 1436355;
    #             """)
    # results = db.cursor.fetchall()
    # for res in results:
    #     print('total_read_count ', res['total_read_count'], ' mutant_allele_read_count ', res['mutant_allele_read_count'], ' type_GT ', res[type_GT])
        


      



if __name__ == '__main__':
    main()