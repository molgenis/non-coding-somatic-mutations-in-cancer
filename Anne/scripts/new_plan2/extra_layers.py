from Database import Database
import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


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
    

def main():
    config = get_config()
    #
    path_db = config['database']  #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    #
    path_file = config['DNase'] #'D:/Hanze_Groningen/STAGE/lagen/2022-04-13_GRCh37_UCNE.bed'
    name_variant = 'DNase'
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    print(len(df_variant))

    add_value(db, name_variant)

    f = open(f'{config["genes_eQTL_etc"]}{name_variant}_num_snps.tsv', 'w') #D:/Hanze_Groningen/STAGE/lagen/
    f.write(f"#Chromosome\tStart\tEnd\tnum_snps_region\n")
    for index, row in df_variant.iterrows():
        num_snps_region = set_value(db, row, name_variant)
        print(index)
        f.write(f"{str(row['#Chromosome'].replace('chr', ''))}\t{int(row['Start'])}\t{int(row['End'])}\t{num_snps_region}\n")
    f.close()
    # Committing the current transactions
    db.mydb_connection.commit()
    # Close database connection
    db.close()

if __name__ == '__main__':
    main()