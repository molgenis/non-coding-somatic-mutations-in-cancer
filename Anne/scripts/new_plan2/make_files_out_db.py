from Database import Database
import sys
import multiprocessing as mp
import pandas as pd

def get_project_ids(db):
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    dict_project = dict()
    for proj in projects:
        dict_project[proj['ID']] = proj['project_ID']
    return dict_project


def test(db, dict_project):
    #https://mungingdata.com/sqlite/create-database-load-csv-python/
    print('BEZIG')
    df = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_project_ID, sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, snp.chr, snp.pos_start, snp.pos_end 
                        FROM project, sum_dosage_GT, snp 
                        WHERE sum_dosage_GT.snp_ID=snp.ID AND sum_dosage_GT.donor_project_ID = project.ID AND (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2) AND sum_dosage_GT.total_read_count_sum >= 33;''', db.mydb_connection)
    # df = pd.read_sql('''SELECT project.cancer, sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, sum_dosage_GT.dosages, snp.in_transcript, snp.in_coding, snp.in_exon, snp.UCNE, snp.TFBS, snp.DNAse, snp.chr, snp.pos_start, snp.pos_end, snp.ref, snp.alt FROM project, sum_dosage_GT, snp WHERE sum_dosage_GT.snp_ID=snp.ID AND sum_dosage_GT.donor_project_ID = project.ID AND (sum_dosage_GT.GT2 = 1 OR sum_dosage_GT.GT2 = 2);''', db.mydb_connection)
    print('HEAD')
    print(df.head())
    print(len(set(df['cancer'])))
    print(set(df['cancer']))
    print(len(set(df['donor_project_ID'])))
    print(set(df['donor_project_ID']))
    print(len(set(df['donor_ID'])))
    print(len(set(df['snp_ID'])))
    print(len(df))
    print(set(df['chr']))
    df = df.loc[df['chr'] != 'MT']
    df['chr'] ='chr' + df['chr'].astype(str)
    print(set(df['chr']))
    df.to_csv('D:/Hanze_Groningen/STAGE/R/all_inf/000ALL.tsv', sep='\t', encoding='utf-8', index=False)

    df1 = df[['chr', 'pos_start', 'pos_end']]
    df1.sort_values(by=['chr', 'pos_start'], inplace=True)
    df1.to_csv('D:/Hanze_Groningen/STAGE/R/select_inf/000ALL.tsv', sep='\t', encoding='utf-8', header=None)

    #breast vs rest
    breast_cancer = df.loc[df['cancer'] == 'Breast']
    breast_cancer.to_csv(f"D:/Hanze_Groningen/STAGE/R/all_inf/breast.tsv", sep='\t', encoding='utf-8', index=False)
    breast_cancer1 = breast_cancer[['chr', 'pos_start', 'pos_end']]
    breast_cancer1.sort_values(by=['chr', 'pos_start'], inplace=True)
    breast_cancer1.to_csv(f"D:/Hanze_Groningen/STAGE/R/select_inf/vs/breast.tsv", sep='\t', encoding='utf-8', header=None)
    #breast vs rest
    nonbreast_cancer = df.loc[df['cancer'] != 'Breast']
    nonbreast_cancer.to_csv(f"D:/Hanze_Groningen/STAGE/R/all_inf/nonbreast.tsv", sep='\t', encoding='utf-8', index=False)
    nonbreast_cancer1 = nonbreast_cancer[['chr', 'pos_start', 'pos_end']]
    nonbreast_cancer1.sort_values(by=['chr', 'pos_start'], inplace=True)
    nonbreast_cancer1.to_csv(f"D:/Hanze_Groningen/STAGE/R/select_inf/vs/nonbreast.tsv", sep='\t', encoding='utf-8', header=None)

    print(nonbreast_cancer1.head())

    # for can in list(set(df['cancer'])):
    #     select_can = df.loc[df['cancer'] == can]
    #     select_can.to_csv(f"D:/Hanze_Groningen/STAGE/R/all_inf/{can.replace(' ', '_')}.tsv", sep='\t', encoding='utf-8', index=False)
    #     select_can1 = select_can[['chr', 'pos_start', 'pos_end']]
    #     select_can1.sort_values(by=['chr', 'pos_start'], inplace=True)
    #     select_can1.to_csv(f"D:/Hanze_Groningen/STAGE/R/select_inf/cancer/{can.replace(' ', '_')}.tsv", sep='\t', encoding='utf-8', header=None)

    

    # for proj in list(set(df['donor_project_ID'])):
    #     select_proj = df.loc[df['donor_project_ID'] == proj]
    #     select_proj.to_csv(f"D:/Hanze_Groningen/STAGE/R/all_inf/{proj}_{dict_project[proj]}.tsv", sep='\t', encoding='utf-8', index=False)
    #     select_proj1 = select_proj[['chr', 'pos_start', 'pos_end']]
    #     select_proj1.sort_values(by=['chr', 'pos_start'], inplace=True)
    #     select_proj1.to_csv(f"D:/Hanze_Groningen/STAGE/R/select_inf/project/{proj}_{dict_project[proj]}.tsv", sep='\t', encoding='utf-8', header=None)


    print('DONE')
    """
    GT = 0 1 2
    15
    40
    2323
    17414723
    18427353
    {'11', '22', '13', '8', '19', '5', '7', '1', '14', '2', '10', '9', '6', 'X', '18', '17', '20', '12', '3', 'Y', '21', '15', '16', 'MT', '4'}
    {'chr6', 'chr8', 'chr15', 'chr11', 'chr10', 'chr7', 'chr20', 'chr2', 'chr3', 'chr14', 'chr4', 'chr1', 'chr12', 'chr18', 'chrX', 'chr5', 'chr16', 'chr19', 'chrY', 'chr17', 'chr22', 'chr21', 'chr13', 'chr9'}
    """
    """
    GT = 1 2
    15
    40
    2271
    7866646
    8126907
    {'19', '6', '8', '1', '3', '17', '21', '16', '10', '22', '15', '4', '20', '14', '7', '2', 'MT', '13', '18', '5', '12', '9', 'X', '11', 'Y'}
    {'chr8', 'chr19', 'chr4', 'chr12', 'chr15', 'chr14', 'chr20', 'chr10', 'chr11', 'chr22', 'chr3', 'chr16', 'chr2', 'chr7', 'chr21', 'chr17', 'chr13', 'chrX', 'chr1', 'chr6', 'chr5', 'chr18', 'chrY', 'chr9'}
    """
           
        #sqlite3 -header -csv /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/db_laatste_copy.db "SELECT sum_dosage_GT.donor_project_ID, sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, snp.chr, snp.pos_start, snp.pos_end FROM sum_dosage_GT, snp WHERE sum_dosage_GT.donor_project_ID = {key} AND sum_dosage_GT.snp_ID=snp.ID;" > /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/test.csv
        # sqlite3 -header -csv /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/db_laatste_copy.db "SELECT sum_dosage_GT.donor_project_ID, sum_dosage_GT.donor_ID, sum_dosage_GT.snp_ID, snp.chr, snp.pos_start, snp.pos_end FROM sum_dosage_GT, snp WHERE sum_dosage_GT.snp_ID=snp.ID;" > /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/test.csv


def main():
    #
    path_db = 'D:/Hanze_Groningen/STAGE/db_laatste_copy.db' #'D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db'
    # Database connection
    db = Database(path_db)
    dict_project = get_project_ids(db)
    print(len(dict_project))
    test(db, dict_project)
    




if __name__ == '__main__':
    main()