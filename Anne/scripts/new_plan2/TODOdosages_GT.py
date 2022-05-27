from Database import Database
import sys
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def set_dosages(db):
    """
    Add a new value to the database: dosages.
    Set dosages to a new value.
    :param db:  The database object
    :return:
    """
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE donor_has_snp
                    ADD `dosages` FLOAT NULL DEFAULT NULL
                    """)
    # set dosages to the correct value
    db.cursor.execute(
            """UPDATE donor_has_snp 
                SET dosages = (CAST(mutant_allele_read_count AS REAL) / CAST(total_read_count AS REAL))
                WHERE total_read_count > 0;""")
    # Add to database
    db.mydb_connection.commit()
    

def set_GT(db, table):
    """
    Add and set GT (genotype) to homozygous alt (2), heterozygous (1) or homozygous ref (0).
    GT remains Nan (or none) when there are no dosages.
    :param db:  The database object
    :return:
    """
    # add dosages
    db.cursor.execute("""
                    ALTER TABLE %s
                    ADD `GT` INT NULL DEFAULT NULL
                    """% table)
    # set GT to the correct value
    # Homozygoot alt
    db.cursor.execute(
            """UPDATE %s 
                SET GT = 2
                WHERE dosages >= 0.8;"""% table)
    # Heterozygoot 
    db.cursor.execute(
            """UPDATE %s 
                SET GT = 1
                WHERE dosages < 0.8 AND dosages > 0.2;"""% table)
    # Homozygoot ref
    db.cursor.execute(
            """UPDATE %s 
                SET GT = 0
                WHERE dosages <= 0.2;"""% table)
    # Add to database
    db.mydb_connection.commit()
    

def set_GT2(db, table):
    """
    Add and set GT2 (genotype).
    The Nan (or none) values are all set to 0 (homozygous ref).
    :param db:  The database object
    :return:
    """
    # add dosages
    db.cursor.execute("""
                    ALTER TABLE %s
                    ADD `GT2` INT NULL DEFAULT NULL
                    """% table)
    # set GT2 to the correct value
    # Make all values of GT2 equal to GT
    db.cursor.execute(
            """UPDATE %s 
                SET GT2 = GT;"""% table)
    # Set all Nan (or none) values to 0 (homozygous ref)
    db.cursor.execute(
            """UPDATE %s 
                SET GT2 = 0
                WHERE GT IS NULL;"""% table)
    db.mydb_connection.commit()

def sum_dosage_GT(db):
    """
    Add new table based on 'donor_had_snp'.
    .....
    :param db:  The database object
    :return:
    """
    # Create new table
    db.cursor.execute("""
        CREATE TABLE sum_dosage_GT AS
            SELECT donor_has_snp.snp_ID, donor_has_snp.donor_ID, donor_has_snp.donor_project_ID,
                    (CAST(SUM(donor_has_snp.mutant_allele_read_count) AS REAL) / CAST(SUM(donor_has_snp.total_read_count) AS REAL)) AS "dosages",
                    SUM(donor_has_snp.total_read_count) AS "total_read_count_sum", 
                    SUM(donor_has_snp.mutant_allele_read_count) AS "mutant_allele_read_count_sum",
                    COUNT(donor_has_snp.total_read_count) AS "number_snps"
            FROM donor_has_snp
            WHERE donor_has_snp.total_read_count > 0
            GROUP BY donor_has_snp.donor_ID, donor_has_snp.snp_ID;
    """)
    # https://database.guide/add-a-foreign-key-to-an-existing-table-in-sqlite/
    db.cursor.execute("""
        PRAGMA foreign_keys = OFF;
    """)
    db.cursor.execute("""
        BEGIN TRANSACTION;
    """)
    # Create new table with forend keys
    db.cursor.execute("""
        CREATE TABLE sum_dosage_GT_new( 
                    `snp_ID` INT NOT NULL,
                    `donor_ID` INT NOT NULL,
                    `donor_project_ID` INT NOT NULL,
                    `dosages` INT NULL DEFAULT NULL,
                    `total_read_count_sum` INT NULL DEFAULT NULL,
                    `mutant_allele_read_count_sum` INT NULL DEFAULT NULL,
                    `number_snps` INT NULL DEFAULT NULL,
                    
                    CONSTRAINT `fk_sum_dosage_GT_new`
                        FOREIGN KEY (`donor_ID`, `donor_project_ID`, `snp_ID`)
                        REFERENCES `donor_has_snp` (`donor_ID`, `donor_project_ID`, `snp_ID`)
                );
    """)
    db.cursor.execute("""
        INSERT INTO sum_dosage_GT_new SELECT * FROM sum_dosage_GT;
    """)
    db.cursor.execute("""
        DROP TABLE sum_dosage_GT;
    """)
    db.cursor.execute("""
        ALTER TABLE sum_dosage_GT_new RENAME TO sum_dosage_GT;
    """)
    db.cursor.execute("""
        COMMIT;
    """)
    db.cursor.execute("""
        PRAGMA foreign_keys = ON;
    """)
    # Committing the current transactions
    db.mydb_connection.commit()


def make_dist_plot_dosages(db, max_donor_id, config):
    # All snps, all donors
    print('All snps, all donors')
    db.cursor.execute("""
                    SELECT dosages, donor_ID, snp_ID, total_read_count_sum, mutant_allele_read_count_sum
                    FROM 'sum_dosage_GT'
                    WHERE total_read_count_sum >= 0 AND mutant_allele_read_count_sum >= 0;
                    """)
    results = db.cursor.fetchall()
    dosages_list = list()
    for res in results:
        dosages_list.append(res['dosages'])
        if res['dosages'] > 1:
            print(res['dosages'], res['donor_ID'], res['snp_ID'], res['total_read_count_sum'], res['mutant_allele_read_count_sum'])
            db.cursor.execute(f"""
                    SELECT dosages, total_read_count, mutant_allele_read_count
                    FROM 'donor_has_snp'
                    WHERE donor_ID = {res['donor_ID']} AND snp_ID = {res['snp_ID']};
                    """)
            check = db.cursor.fetchall()
            for ch in check:
                print('--------------------', ch['total_read_count'], ch['mutant_allele_read_count'])
    # Using filter() method to filter None values
    filtered_list = list(filter(None, dosages_list))
    print(max(filtered_list), min(filtered_list))
    sns.displot(dosages_list)
    plt.tight_layout()
    plt.savefig(f"{config['genes_eQTL_etc']}dosages_snps_donor.png")
    plt.clf()
    plt.close()
    # # All donors
    print('All donors')
    db.cursor.execute("""
                    SELECT donor_ID
                    FROM 'sum_dosage_GT'
                    WHERE total_read_count_sum >= 0 and mutant_allele_read_count_sum >= 0;
                    """)
    results = db.cursor.fetchall()
    # dosages_dict = dict()
    count_snps = dict()
    for res in results:
        if res['donor_ID'] in count_snps:
            # dosages_dict[res['donor_ID']].append(res['dosages'])
            count_snps[res['donor_ID']] += 1
        else:
            # dosages_dict[res['donor_ID']] = [res['dosages']]
            count_snps[res['donor_ID']] = 1
    print('make plots donors')

    count_snps = dict(sorted(count_snps.items()))
    plt.bar(count_snps.keys(), count_snps.values(), color='g')
    plt.tight_layout()
    plt.savefig(f"{config['genes_eQTL_etc']}count_snps.png")
    plt.clf()
    plt.close()
    # for key, value in dosages_dict.items():        
    #     sns.displot(value)
    #     plt.tight_layout()
    #     plt.savefig(f"D:/Hanze_Groningen/STAGE/eQTL/dosages_{key}donor.png")
    #     plt.clf()
    #     plt.close()
    db.cursor.execute("""
                    SELECT donor_ID
                    FROM 'sum_dosage_GT';
                    """)
    results = db.cursor.fetchall()
    # dosages_dict2 = dict()
    list_ID = list()
    count_snps2 = dict()
    for res in results:
        if res['donor_ID'] in count_snps2:
            # dosages_dict2[res['donor_ID']].append(res['dosages'])
            count_snps2[res['donor_ID']] += 1
        else:
            # dosages_dict2[res['donor_ID']] = [res['dosages']]
            count_snps2[res['donor_ID']] = 1
            list_ID.append(res['donor_ID'])

    

    count_snps2 = dict(sorted(count_snps2.items()))
    plt.bar(count_snps2.keys(), count_snps2.values(), color='g')
    plt.tight_layout()
    plt.savefig(f"{config['genes_eQTL_etc']}count_snps2_zonderfilter.png")
    plt.clf()
    plt.close()

    list_ID_sort = list(filter(None, list_ID))    
    range_1 = list(range(1, max_donor_id+1))

    listtest = list(set(range_1) - set(list_ID_sort))
    print(listtest)
    print(len(listtest))
    return listtest
    
        
def make_dist_plot_tot(db, config):
    # All snps, all donors
    print('All snps, all donors')
    db.cursor.execute("""
                    SELECT total_read_count_sum
                    FROM 'sum_dosage_GT'
                    WHERE total_read_count_sum >= 0 AND mutant_allele_read_count_sum >= 0;
                    """)
    results = db.cursor.fetchall()
    tot_list = list()
    for res in results:
        tot_list.append(res['dosages'])
    sns.displot(tot_list)
    plt.tight_layout()
    plt.savefig(f"{config['genes_eQTL_etc']}tot_snps_donor.png")
    plt.clf()
    plt.close()


def get_min_max_snpID(db):
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
                    FROM donor;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        max_donor_id = res[0]
    # MIN
    db.cursor.execute("""
                    SELECT MIN(ID) 
                    FROM donor;
                    """)
    results = db.cursor.fetchall()
    for res in results:
        min_donor_id = res[0]
    return max_donor_id, min_donor_id




def main():
    config = get_config()
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" 
    # Database connection
    db = Database(path_db)
    # set_dosages(db)
    set_GT(db, 'donor_has_snp')
    set_GT2(db, 'donor_has_snp')
    # sum_dosage_GT(db)
    set_GT(db, 'sum_dosage_GT')
    set_GT2(db, 'sum_dosage_GT')

    max_donor_id, min_donor_id = get_min_max_snpID(db)

    # make_dist_plot_tot(db, config)
    myfile = open(config['Number_and_none_values'], 'w')
    myfile.writelines(f'donor\tcountNan\tcountNum\tcountNan_sum\tcountNum_sum\n')

    for value in range(1, max_donor_id+1): #listtest:
        countNan = 0
        countNum = 0
        # print(f'------value: {value}')        
        db.cursor.execute(f"""
                        SELECT snp_ID, total_read_count, mutant_allele_read_count
                        FROM donor_has_snp
                        WHERE donor_ID = {value};
                        """)
        results = db.cursor.fetchall()
        for res in results:
            if res['total_read_count'] is None:
                countNan += 1
            else:
                countNum += 1
        countNan_sum = 0
        countNum_sum = 0
        db.cursor.execute(f"""
                        SELECT snp_ID, total_read_count_sum, mutant_allele_read_count_sum
                        FROM sum_dosage_GT
                        WHERE donor_ID = {value};
                        """)
        results = db.cursor.fetchall()
        for res in results:
            if res['total_read_count_sum'] is None:
                countNan_sum += 1
            else:
                countNum_sum += 1
        myfile.writelines(f'{value}\t{countNan}\t{countNum}\t{countNan_sum}\t{countNum_sum}\n')

    myfile.close()
    
      



if __name__ == '__main__':
    main()