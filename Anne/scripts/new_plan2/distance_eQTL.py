import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def search_snp(db, row, region):
    """
    
    """
    # Count snps in region
    db.cursor.execute(
        """SELECT ID, pos_start, chr
            FROM snp
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(row['SNPChr']), int(row['SNPPos']-region), int(row['SNPPos']+region)))
    
    results = db.cursor.fetchall()
    return results

def closest_snp(results, row):
    """
    
    """
    for res in results:
        if res['pos_start'] < distance_eQTL:
            distance_eQTL = res['pos_start']
            distance_snp_ID = f"{row['ID']}_{row['chr']}_{row['pos_start']}"
    return distance_snp_ID, distance_eQTL

def change_range(db, row, region):
    """
    
    """
    results = search_snp(db, row, region)
    if len(results) > 0:
        distance_snp_ID, distance_eQTL = closest_snp(results, row)
        return distance_snp_ID, distance_eQTL
    else:
        return 0, ''


def search_close_snp(db, df_strong_eQTL, config):
    """
    
    """
    f = open(config['distance_eqtl_path'], 'w') #D:/Hanze_Groningen/STAGE/eQTL/
    f.write(f"SNP_eQTL\chr\pos\Gene\tsnp_ID\tchr_snp\tpos_snp\tdistance\n")
    for index, row in df_strong_eQTL.iterrows():
        distance_eQTL = 100000
        for region in range(1000, 1000000, 1000):
            distance_snp_ID, distance_eQTL = change_range(db, row, region)
            if distance_snp_ID != 0:
                break
        if distance_snp_ID != 0:
            info_distance_snp_ID = distance_snp_ID.split('_')
            f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t{info_distance_snp_ID[0]}\t{info_distance_snp_ID[1]}\t{info_distance_snp_ID[2]}\t{distance_eQTL}\n")
        else:
            f.write(f"{str(row['SNP'])}\t{int(row['SNPChr'])}\t{int(row['SNPPos'])}\t{row['Gene']}\t-\t-\t-\t-\n")
    f.close()


    # results = search_snp(db, row, 100)
        # if len(results) > 0:
        #     distance_snp_ID, distance_eQTL = closest_snp(results, row)
        # else:
        #     results = search_snp(db, row, 500)
        #     if len(results) > 0:
        #         distance_snp_ID, distance_eQTL = closest_snp(results, row)
        #     else:
        #         results = search_snp(db, row, 1000)
        #         if len(results) > 0:
        #             distance_snp_ID, distance_eQTL = closest_snp(results, row)
        #         else:
        #             pass  

def main():
    config = get_config('gearshift')
    #
    path_db = config['database'] #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    path_strong_eQTL = config['strong_eqtl_path'] #"D:/Hanze_Groningen/STAGE/eQTL/eqtl_v1013_lead_snp_gene_with_info.txt"
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t')
    search_close_snp(db, df_strong_eQTL, config)




if __name__ == '__main__':
    main()


# loop over file met sterkste eQTL
# zoek 100 bp voor en na een eQTL
# snps
    # ja
        # loop over snps en kijk welke de kleinste afstand heeft tot eQTL
    # nee
        # Maak van de afstand 500
        # etc