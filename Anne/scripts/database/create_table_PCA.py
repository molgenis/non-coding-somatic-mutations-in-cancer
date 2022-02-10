import pandas as pd
from collections import defaultdict
from Database import Database
from multiprocessing import Pool
import multiprocessing as mp
import time
import sys




def filter_donor(test_file, dict_test, donor, chr, i_start, i): 
    filter_file = test_file[(test_file['donor_id'] == str(donor)) & (test_file['chr'] == str(chr)) & 
                                        (test_file['pos_start'] >= int(i_start+1)) & (test_file['pos_end'] <= int(i)) &
                                        (test_file['seq_strategy'] == 'WGS')]
    if len(filter_file) != 0:
        if donor in dict_test:
            dict_test[donor].append(len(filter_file))
        else:
            dict_test[donor] = [len(filter_file)]
        print(donor)
        print(len(filter_file))
    else:
        if donor in dict_test:
            dict_test[donor].append(int(0))
        else:
            dict_test[donor] = [int(0)]
    return dict_test

def make_table(donors, test1_file, chr, chr_length, steps, dict_donor_project, path_save):
    print(chr)
    dict_test = dict()
    test_file = test1_file[test1_file['chr'] == str(chr)]
    i_start=0
    for i in range(0,chr_length[chr],steps):
        print(f'chr{chr}:{i_start+1} - {i}')
        if i != 0:
            if 'index' in dict_test:
                add_value = f'chr{chr}:{i_start+1}-{i}'
                dict_test['index'].append(add_value)
            else:
                dict_test['index'] = [f'chr{chr}:{i_start+1}-{i}']
            # cursor.execute("""SELECT *
            #                 FROM donor""")
            # donors = cursor.fetchall()
            for donor in donors:
                dict_test = filter_donor(test_file, dict_test, donor, chr, i_start, i)
                if chr_length[chr] == (i + (chr_length[chr] % steps)):
                    print(i)
                    last_i = (i + (chr_length[chr] % steps))
                    print(last_i)
                    print(chr_length[chr])
                    dict_test = filter_donor(test_file, dict_test, donor, chr, i_start, last_i)
        i_start = i
    
    df_chr = pd.DataFrame.from_dict(dict_test,orient='index').transpose() #pd.DataFrame([dict_test])
    testje = pd.DataFrame.from_dict(dict_donor_project,orient='index').transpose()
    result = pd.concat([testje, df_chr])
    result.to_csv(f'{path_save}chr{chr}.tsv', sep='\t', encoding='utf-8', index=False) 


def main():
    start_time = time.time()
    steps= 1000000
    #https://en.wikipedia.org/wiki/Human_genome
    chr_length = {'17':83257441, '18':80373285, '19':58617616, '20':64444167} #{'1':248956422, '2':242193529, '3':198295559, '4':190214555,
                # '5':181538259, '6':170805979, '7':159345973, '8':145138636,
                # '9':138394717, '10':133797422, '11':135086622, '12':133275309,
                # '13':114364328, '14':107043718, '15':101991189, '16':90338345,
                # '17':83257441, '18':80373285, '19':58617616, '20':64444167,
                # '21':46709983, '22':50818468, 'X':156040895, 'Y':57227415} #{'17':83257441, '18':80373285, '19':58617616, '20':64444167}
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db" #sys.argv[1]
    path_save = "D:/Hanze_Groningen/STAGE/UMAP/" #sys.argv[2]
    # Database connection
    db = Database(path_db)
    #project
    db = Database(path_db)
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    dict_project = dict()
    for proj in projects:
        dict_project[proj['ID']] = proj['project_ID']
    #donor
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    donor_list =  []
    dict_donor_project = {'index':'project_ID'}
    for donor in donors:
        donor_list.append(donor['donor_ID'])
        dict_donor_project[donor['donor_ID']] = dict_project[donor['project_ID']]

    


    file_path = "D:/Hanze_Groningen/STAGE/db/files/ALL-US_db.tsv" #sys.argv[3]
    test1_file = pd.read_csv(file_path, sep='\t')
    print(test1_file.dtypes)

    print(mp.cpu_count()-1)
    arg_multi_list = []
    for key, value in chr_length.items():
        arg_multi_list.append((donor_list, test1_file, key, chr_length, steps, dict_donor_project, path_save))

    pool = Pool(processes=mp.cpu_count()-1)
    pool.starmap(func=make_table, iterable=arg_multi_list)
    pool.close()
    pool.join()



    # make_table(donors, test1_file, '1', chr_length, steps)
    db.close()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    main()
