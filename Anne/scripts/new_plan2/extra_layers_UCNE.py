from Database import Database
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from collections import Counter



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
    # db.cursor.execute(
    #     """UPDATE snp
    #         SET %s = TRUE
    #         WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
    #     (str(name_variant), str(row['#Chromosome'].replace('chr', '')), int(row['Start']), int(row['End'])))
    # Count snps in region
    db.cursor.execute(
        """SELECT ID
            FROM snp
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(row['#Chromosome'].replace('chr', '')), int(row['Start']), int(row['End'])))
    
    results = db.cursor.fetchall()
    return len(results)

def close_to(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, donor_cancer_list, filter_num):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (gene_file).
    :param db:                The database object
    :param gene:              The name of the gene
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param gene_file:         The file after which the genes with their specific region before or after
                              a gene are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param gene_name_list --> gene_region_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param sparse_matrix_region:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return: sparse_matrix_region:   A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
             gene_name_list:  List of gene names (to be used later as columns in the sparse matrix)
    """
    # # Get index of the gene name in the list
    # gene_index = gene_name_list.index(gene)
    # if with_type == 'genes':
    #     # Replace the name with gene_name:chromosoom:start_pos-end_pos
    #     gene_name_list[gene_index] = f'{gene}:chr{chr}:{int(start_pos)}-{int(end_pos)}'
    
    #
    donor_read_count = dict()
    # Make donor_list_snp. List because double donors count
    donor_list_snp = list()
    #
    snp_ID_list = list()
    
    
    db.cursor.execute("""
                    SELECT sum_dosage_GT.snp_ID , sum_dosage_GT.donor_ID, 
                           sum_dosage_GT.total_read_count_sum , sum_dosage_GT.mutant_allele_read_count_sum,
                           sum_dosage_GT.dosages
                    FROM snp, sum_dosage_GT
                    WHERE snp.chr = '%s' AND snp.pos_start >= %s AND snp.pos_end <= %s AND sum_dosage_GT.total_read_count_sum > %s 
                            AND sum_dosage_GT.mutant_allele_read_count_sum > 0 AND snp.ID = sum_dosage_GT.snp_ID
                    GROUP BY sum_dosage_GT.snp_ID, sum_dosage_GT.donor_ID;
                    """ %
                    (str(chr), int(start_pos), int(end_pos), int(filter_num)))

    results = db.cursor.fetchall()
    if len(results) > 0:
        for res in results:
            print(f'{res[0]} - {res[1]} - {res[2]} - {res[3]} - {res[4]}')
            donor_list_snp.append(donor_dict[res['donor_ID']])
            snp_ID_list.append(res['snp_ID'])
            # donor_index = donor_list.index(res['sum_dosage_GT.donor_ID'])
            # sparse_matrix_overall[donor_index, gene_index] += 1
            if donor_dict[res['donor_ID']] in donor_read_count:
                donor_read_count[donor_dict[res['donor_ID']]][0] += res['total_read_count_sum']
                donor_read_count[donor_dict[res['donor_ID']]][1] += res['dosages']
                donor_read_count[donor_dict[res['donor_ID']]][2] += 1
            else:
                donor_read_count[donor_dict[res['donor_ID']]] = [res['total_read_count_sum'], res['dosages'], 1]
            # Loop over dict: donor_read_count
        for key, value in donor_read_count.items():
            donor_index = donor_list.index(key)
        # Make cancer_list
        cancer_list = list()
        for donor in donor_list_snp:
            donor_index = donor_list.index(donor)
            # Adds +1 to the sparse_matrix at the position of donor region
            # sparse_matrix_overall[donor_index, gene_index] += 1
            cancer_list.append(donor_cancer_list[donor_index])
        # Creates a dictionary from the list with as key the name (cancer type or donor ID) and
        # as value how often that name occurs in the list.
        cancer_count = dict(Counter(cancer_list))
        donor_count = dict(Counter(donor_list_snp))
    
        gene_file.write(str(filter_num) + '\t' + str(chr) + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
        len(snp_ID_list)) + '\t' + ','.join(map(str, snp_ID_list)) + '\t' + str(len(donor_list_snp)) + '\t' + str(
        donor_count) + '\t' + str(cancer_count) + '\n')
    else:       
        gene_file.write(str(filter_num) + '\t' + str(chr) + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
        len(snp_ID_list)) + '\t-\t-\t-\t-\n')
    


def make_file_extra_layers(df_variant, chr, name_variant, path_db, config, donor_dict, donor_list, donor_cancer_list):
    # Database connection
    db = Database(path_db)
    # f = open(f'{config["genes_eQTL_etc"]}{name_variant}_{chr}_num_snps.tsv', 'w') #D:/Hanze_Groningen/STAGE/lagen/
    # f.write(f"#Chromosome\tStart\tEnd\tnum_snps_region\n")
     # The header for the files before_gene_file and after_gene_file
    header_file = 'filter\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    
    # Make file
    f = open(f'{config["genes_eQTL_etc"]}{name_variant}_{chr}_num_snps_ALL.tsv', 'w')
    f.write(header_file)
    for index, row in df_variant.iterrows():
        # num_snps_region = set_value(db, row, name_variant)
        close_to(db, index, row['#Chromosome'].replace('chr', ''), row['Start'], row['End'], f, donor_dict, donor_list, donor_cancer_list, 33)
        print(index)
        # f.write(f"{str(row['#Chromosome'].replace('chr', ''))}\t{int(row['Start'])}\t{int(row['End'])}\t{num_snps_region}\n")
    f.close()
    # Committing the current transactions
    db.mydb_connection.commit()
    

def main():
    config = get_config()
    #
    path_db = config['database_UCNE'] #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    #
    path_file = config['UCNE'] #'D:/Hanze_Groningen/STAGE/lagen/2022-04-13_GRCh37_UCNE.bed'
    name_variant = 'UCNE'
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    print(len(df_variant))

    # add_value(db, name_variant)

    arg_multi_list = []
    for item in list(set(df_variant['#Chromosome'])):
        select_variant_df = df_variant.loc[df_variant['#Chromosome'] == item]
        print(item)
        arg_multi_list.append((select_variant_df, item, name_variant, path_db, config, donor_dict, donor_list, donor_cancer_list))

    pool = Pool(processes=mp.cpu_count())
    pool.starmap(func=make_file_extra_layers, iterable=arg_multi_list)
    pool.close()
    pool.join()

    # Close database connection
    db.close()

if __name__ == '__main__':
    main()