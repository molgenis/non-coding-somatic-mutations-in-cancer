from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix
# from matplotlib import pyplot as plt



def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/db_laatste_copy.db" #"/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydb_L.db"  # /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    print('ZOEk')

    # db.cursor.execute("""
    #                 SELECT *
    #                 FROM snp
    #                 WHERE UCNE = 1
    #                 """)
    # print('KLAAR')
    # results = db.cursor.fetchall()
    # print(len(results))
    # for res in results:
    #     print(res['ID'])
    print('BEZIG')
    df = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                            sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, 
                            snp.chr, snp.pos_start, snp.pos_end , snp.DNase, snp.TFBS, snp.UCNE
                    FROM project, sum_dosage_GT, snp 
                    WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                              sum_dosage_GT.donor_project_ID = project.ID AND 
                              (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                              sum_dosage_GT.total_read_count_sum >= 33 AND snp.UCNE = 1;''', db.mydb_connection)
    print('KLAAR')
    print(df.head())
    print(set(df['TFBS']))
    print(set(df['UCNE']))
    print(set(df['DNase']))


                              
    
    # # Path of the genes and there positions
    # gene_path = 'D:/Hanze_Groningen/STAGE/db/all_genes_new.tsv'#"/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/all_genes_new.tsv"  # snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed

    # # Read gene file
    # gene_df = pd.read_csv(gene_path, sep='\t')

    # # Region before the start position of a gene or after the stop position of a gene
    # position_out_gene = 2000
    # # Region after the start position of a gene or before the stop position of a gene
    # position_in_gene = 250

    # for index, row in gene_df.iterrows():
    #     chr = row['hg19.knownGene.chrom'].replace('chr', '')
    #     print('##################################')
    #     print(f"{row['hg19.knownGene.txStart'] - position_out_gene} - {row['hg19.knownGene.txStart'] + position_in_gene}")
        
    #     db.cursor.execute("""
    #                     SELECT snp.ID, sum_dosage_GT.snp_ID , sum_dosage_GT.donor_ID, 
    #                            sum_dosage_GT.total_read_count_sum , sum_dosage_GT.mutant_allele_read_count_sum, 
    #                            snp.pos_start, snp.pos_end
    #                     FROM snp, sum_dosage_GT
    #                     WHERE snp.chr = '%s' AND snp.pos_start >= %s AND snp.pos_end <= %s AND sum_dosage_GT.total_read_count_sum > 0 
    #                           AND sum_dosage_GT.mutant_allele_read_count_sum > 0 AND snp.ID = sum_dosage_GT.snp_ID
    #                     GROUP BY sum_dosage_GT.snp_ID, sum_dosage_GT.donor_ID;
    #                     """ %
    #                     (str(chr), int(row['hg19.knownGene.txStart'] - position_out_gene), int(row['hg19.knownGene.txStart'] + position_in_gene)))

    #     results = db.cursor.fetchall()
    #     for res in results:
    #         print(f'{res[0]} - {res[1]} - {res[2]} - {res[3]} - {res[4]} - {res[5]} - {res[6]}')

        # db.cursor.execute("""
        #                 SELECT ID
        #                 FROM snp
        #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;
        #                 """ %
        #                 (str(chr), int(row['hg19.knownGene.txStart'] - position_out_gene), int(row['hg19.knownGene.txStart'] + position_in_gene)))

        # results = db.cursor.fetchall()
        # for res in results:
        #     db.cursor.execute("""
        #                 SELECT snp_ID, donor_ID, total_read_count, dosages
        #                 FROM 'donor_has_snp'
        #                 WHERE snp_ID = %s AND total_read_count > %s AND mutant_allele_read_count > 0;
        #                 """ % (res[0], 0))
        #     donors = db.cursor.fetchall()
        #     for don in donors:
        #         print(f'{res[0]} - {don[0]} - {don[1]} - - {don[2]}')



if __name__ == '__main__':
    main()

# def test():
#     count = 0
#     for i in range(0,10):
#         print(i)
#         for x in range(1,4):
#             print('xx', x)
#             if ( i%2 ) == 0 and (x%2) == 0:
#                 print('jaaaa')
#                 break
    
#     print(count)

# test()
    




