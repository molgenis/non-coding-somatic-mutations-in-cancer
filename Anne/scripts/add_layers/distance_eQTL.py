#!/usr/bin/env python3

# Imports
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config



def multiprocess_search_close_snp(df_strong_eQTL, region, path_save_file, config):
    """
    Sees which SNPs are within a specific region of an eQTL
    :param df_strong_eQTL: Df with strongest eQTLs
    :param region:         Region from which it looks at a distance from an eQTL if there are SNPs within that region
    :param path_save_file: Path where the file is saved
    :param config:         Dictionary with as keys the name of the paths and as value the paths           
    :return: 
    """
    path_db = config['database']
    # Database connection
    db = Database(path_db)
    f = open(path_save_file, 'w')
    f.write(f"SNP_eQTL\tchr\tpos\tGene\tsnp_ID\tchr_snp\tpos_snp\tdistance\n")
    # Loop over file
    for index, row in df_strong_eQTL.iterrows():
        db.cursor.execute(
        """SELECT ID, chr, pos_start, abs(CAST(%s AS REAL) - CAST(pos_start AS REAL)) AS distance_eQTL
            FROM snp
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
            ORDER BY distance_eQTL
            LIMIT 1;""" %
        (int(row['SNPPos']), str(row['SNPChr']), int(row['SNPPos']-region), int(row['SNPPos']+region)))
    
        results = db.cursor.fetchall()
        if len(results) != 0:
            for res in results:
                f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}\n")
        else:
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
    # Call get_config
    config = get_config('gearshift')
    # Get path
    path_strong_eQTL = config['strong_eqtl_path']
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t')
    # Region from which it looks at a distance from an eQTL if there are SNPs within that region
    region = 10000
    
    cpus = mp.cpu_count()
    list_paths_save = list()
    arg_multi_list = []
    for i in range(0, len(df_strong_eQTL)+1000, 1000):
        df = df_strong_eQTL.iloc[i:i+1000]
        path_save_file = f'{config["genes_eQTL_etc"]}distance_eqtl_{i}.tsv'
        list_paths_save.append(path_save_file)
        # add parameters needed to run multiprocess_search_close_snp function
        arg_multi_list.append((df, region, path_save_file, config))
    # Multiprocess
    pool = Pool(processes=cpus)
    pool.starmap(func=multiprocess_search_close_snp, iterable=arg_multi_list)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()