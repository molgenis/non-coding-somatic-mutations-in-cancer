import pandas as pd
from collections import defaultdict
from Database import Database
from multiprocessing import Pool
import multiprocessing as mp
import time




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

def make_table(donors, test1_file, chr, chr_length, steps):
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
    df_chr.to_csv(f'D:/Hanze_Groningen/STAGE/UMAP/chr{chr}.tsv', sep='\t', encoding='utf-8')


def main():
    start_time = time.time()
    steps= 1000000
    #https://en.wikipedia.org/wiki/Human_genome
    chr_length = {'1':248956422, '2':242193529, '3':198295559, '4':190214555,
                '5':181538259, '6':170805979, '7':159345973, '8':145138636,
                '9':138394717, '10':133797422, '11':135086622, '12':133275309,
                '13':114364328, '14':107043718, '15':101991189, '16':90338345,
                '17':83257441, '18':80373285, '19':58617616, '20':64444167,
                '21':46709983, '22':50818468, 'X':156040895, 'Y':57227415} #{'17':83257441, '18':80373285, '19':58617616, '20':64444167}
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db"
    # Database connection
    db = Database(path_db)
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    donor_list =  []
    for donor in donors:
        donor_list.append(donor['donor_ID'])


    file_path = "D:/Hanze_Groningen/STAGE/db/files/ALL-US_db.tsv"
    test1_file = pd.read_csv(file_path, sep='\t')
    print(test1_file.dtypes)

    print(mp.cpu_count()-1)
    # arg_multi_list = []
    # for key, value in chr_length.items():
    #     arg_multi_list.append((donor_list, test1_file, key, chr_length, steps))

    # pool = Pool(processes=mp.cpu_count()-1)
    # pool.starmap(func=make_table, iterable=arg_multi_list)
    # pool.close()
    # pool.join()



    make_table(donors, test1_file, '1', chr_length, steps)
    db.close()
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    main()


# dict_test = dict()
# for chr in [1]:
    # test_file = test1_file[test1_file['chr'] == str(chr)]
    # i_start=0
    # for i in range(0,chr_length[chr],steps):
    #     count = 0
    #     print(f'{i_start+1} - {i}')
    #     if i != 0:
    #         if 'index' in dict_test:
    #             add_value = f'chr{chr}:{i_start+1}-{i}'
    #             dict_test['index'].append(add_value)
    #         else:
    #             dict_test['index'] = [f'chr{chr}:{i_start+1}-{i}']


    #         cursor.execute("""SELECT *
    #                         FROM donor""")
    #         donors = cursor.fetchall()
    #         for donor in donors:
    #             dict_test = filter_donor(test_file, dict_test, donor, chr, i_start, i)
    #         if i == (chr_length[chr] + (chr_length[chr] % steps)):
    #             print(i)
    #             last_i = (i + (chr_length[chr] % steps))
    #             print(last_i)
    #             print(chr_length[chr])
    #             dict_test = filter_donor(test_file, dict_test, donor, chr, i_start, last_i)
    #     i_start = i

# DATAF = pd.DataFrame([dict_test])
# print(DATAF)
# print(len(DATAF.columns))
    
                #     print('NOOO')
                # print(dict_test)

        # cursor.execute("""SELECT *
        #                 FROM donor""")
        # donors = cursor.fetchall()
        # for donor in donors:
        #     cursor.execute("""SELECT *
        #                     FROM donor_has_snp
        #                     WHERE donor_ID = %s"""% 
        #                     (int(donor['ID'])))
        #     all_snps = cursor.fetchall()
        #     for snp in all_snps:
        #         cursor.execute("""SELECT COUNT(*)
        #                         FROM snp
        #                         WHERE ID = %s AND seq_strategy = 'WGS' AND chr = %s 
        #                         AND pos_start >= %s AND pos_start <= %s"""% 
        #                         (int(snp['snp_ID']), str(chr), int(i_start+1), int(i)))
        #         select_snp = cursor.fetchall()
        #         for select in select_snp:
        #             if select[0] != 0:
        #                 print(select[0])
            

        # i_start = i
    # if (chr_length[chr] % steps) != 0:
    #     print(chr_length[chr] % steps)
    #     print(f'{i_start +1} - {i + (chr_length[chr] % steps)}')