from Database import Database
import pandas as pd
from multiprocessing import Pool
import multiprocessing as mp
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config
from collections import Counter


def layer_run(df_variant, name_variant, df, path_save, chr):    
    header_file = 'filter\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                '\tdonor_count\tcancer_count\n'    
    # Make file
    f = open(f'{path_save}{name_variant}_chr{chr}_num_snps_ALL_2.tsv', 'w')
    f.write(header_file)
    for index, row in df_variant.iterrows():
        print(index)
        x = df[(df['chr']==row['#Chromosome']) & (df['pos_start']>=row['Start']) & (df['pos_end']<=row['End'])]
        donor_list = list(x['donor_ID'])
        donor_count = dict(Counter(donor_list))
        cancer_list = list(x['cancer'])
        cancer_count = dict(Counter(cancer_list))
        if len(x) > 0:
            f.write(str(33) + '\t' + str(row['#Chromosome']) + '\t' + str(row['Start']) + '\t' + str(row['End']) + '\t' + str(
            len(x['snp_ID'])) + '\t' + ','.join(map(str, list(x['snp_ID']))) + '\t' + str(len(donor_list)) + '\t' + str(
            donor_count) + '\t' + str(cancer_count) + '\n')
        else:
            f.write(str(33) + '\t' + str(row['#Chromosome']) + '\t' + str(row['Start']) + '\t' + str(row['End']) + '\t-\t-\t-\t-\t-\n')

    f.close()


def run_code(config, df, chr):
    # UCNE
    path_file = config['UCNE']
    df_variant = pd.read_csv(path_file, sep='\t', compression='gzip')
    df_variant['#Chromosome'] = df_variant['#Chromosome'].str.replace('chr', '')
    print(len(df_variant))
    name_variant = 'UCNE'
    path_save = config["genes_eQTL_etc"]
    layer_run(df_variant, name_variant, df, path_save, chr)


def main():
    config = get_config()
    path_db = config['database_UCNE']  #"/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydb_L.db"  # /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)

    df = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, 
                            sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, 
                            snp.chr, snp.pos_start, snp.pos_end , snp.DNase, snp.TFBS, snp.UCNE
                    FROM project, sum_dosage_GT, snp 
                    WHERE sum_dosage_GT.snp_ID=snp.ID AND 
                                sum_dosage_GT.donor_project_ID = project.ID AND 
                                (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND 
                                sum_dosage_GT.total_read_count_sum >= 33 ;''', db.mydb_connection)
    chromosomes_list = list(set(df['chr']))
    
    cpus = mp.cpu_count()
    arg_multi_list = []
    for chrom in chromosomes_list:
        df_chr = df[df['chr'] == chrom]
        arg_multi_list.append((config, df_chr, chrom))

    pool = Pool(processes=cpus)
    pool.starmap(func=run_code, iterable=arg_multi_list)
    pool.close()
    pool.join()



if __name__ == '__main__':
    main()
