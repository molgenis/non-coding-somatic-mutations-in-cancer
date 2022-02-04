#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database
#TRUE = 1 and FALSE=0


def close_eQTL(position_gene, cursor, where, chr, start, end, count_number, set_snps, count_unique):    

    # START
    cursor.execute("""
                    SELECT *
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND %s
                    """ %
        (str(chr), int(start)-position_gene, int(start), str(where)))
    results = cursor.fetchall()
    start_set = set()
    for res in results:
        start_set.add(res['ID'])

    count_number += len(list(start_set))
    # letters in a but not in set_snps
    unique_start = list(start_set - set(set_snps))
    count_unique += len(unique_start)
    set_snps.update(unique_start)
    # if count_number > 0:
    #     print('---NEWSTART---')
    #     print(count_number)
    #     print(set_snps)
    #     print(count_unique)


    # END
    cursor.execute("""
                    SELECT *
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND %s
                    """ %
        (str(chr), int(end), int(end)+position_gene, str(where)))
    results = cursor.fetchall()
    end_set = set()
    for res in results:
        end_set.add(res['ID'])

    count_number += len(end_set)
    # letters in a but not in set_snps
    unique_end = list(end_set - set(set_snps))
    count_unique += len(unique_end)
    set_snps.update(unique_end)
    # if count_number > 0:
    #     print('---NEWEND---')
    #     print(count_number)
    #     print(set_snps)
    #     print(count_unique)
    if count_number > 0:
        print(count_number)
        # print(set_snps)
        print(count_unique)
        print('-------')

    return count_number, set_snps, count_unique


#'', 100, db.cursor, 'in_transcript', 'snp', where+f'in_transcript = {value}', 'txStart', 'txEnd'
def place_eQTL(eQTL_path, position_gene, cursor, where):
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO
    eQTL_df = pd.read_csv(eQTL_path, sep='\t')
    count_number = 0
    count_unique = 0
    set_snps = set()
    # CLOSE FROM GENES
    for index, row in eQTL_df.iterrows():
        count_number, set_snps, count_unique = close_eQTL(position_gene, cursor, where, row['SNPChr'], row['SNPPos'], row['SNPPos'], count_number, set_snps, count_unique)
    print('END')
    print(count_number)