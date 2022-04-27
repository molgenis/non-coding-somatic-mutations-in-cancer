from Database import Database
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp


def add_value(db, name_variant):
    """
    Adds value to the database (table snp).
    :param db:  The database object
    :return:
    """
    # Add in_transcript
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD {name_variant} BOOLEAN DEFAULT(FALSE)
                    """)
    
    # Committing the current transactions
    db.mydb_connection.commit()

def set_value(db, row, name_variant):
    """
    
    """
    print(row['#Chromosome'], row['Start'], row['End'])
    # Update in_transcript
    db.cursor.execute(
        """UPDATE snp
            SET %s = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(name_variant), str(row['#Chromosome'].replace('chr', '')), int(row['Start']), int(row['End'])))
    # Count snps in region
    db.cursor.execute(
        """SELECT ID
            FROM snp
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(row['#Chromosome'].replace('chr', '')), int(row['Start']), int(row['End'])))
    
    results = db.cursor.fetchall()
    return len(results)

def make_file_extra_layers(df_variant, chr, name_variant, path_db):
    # Database connection
    db = Database(path_db)
    f = open(f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/genes_eQTL_etc/{name_variant}_{chr}_num_snps.tsv', 'w') #D:/Hanze_Groningen/STAGE/lagen/
    f.write(f"#Chromosome\tStart\tEnd\tnum_snps_region\n")
    for index, row in df_variant.iterrows():
        num_snps_region = set_value(db, row, name_variant)
        print(index)
        f.write(f"{str(row['#Chromosome'].replace('chr', ''))}\t{int(row['Start'])}\t{int(row['End'])}\t{num_snps_region}\n")
    f.close()
    # Committing the current transactions
    db.mydb_connection.commit()
    

def main():
    #
    path_db = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/db_laatste_copy.db' #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    #
    path_file = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/GREEN_DB/2022-04-13_GRCh37_DNase.merged.bed.gz' #'D:/Hanze_Groningen/STAGE/lagen/2022-04-13_GRCh37_UCNE.bed'
    name_variant = 'DNase'
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    print(len(df_variant))

    # add_value(db, name_variant)

    arg_multi_list = []
    for item in list(set(df_variant['#Chromosome'])):
        select_variant_df = df_variant.loc[df_variant['#Chromosome'] == item]
        print(item)
        arg_multi_list.append((select_variant_df, item, name_variant, path_db))

    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=make_file_extra_layers, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # Close database connection
    db.close()

if __name__ == '__main__':
    main()