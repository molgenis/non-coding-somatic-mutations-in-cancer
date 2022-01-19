#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database

def close_eQTL(cursor, position_gene, start, end, chr, where, count_number, set_snps, count_unique):
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


def place_eQTL(position_gene, cursor, where):
    cursor.execute("""
                    SELECT *
                    FROM 'snp'
                    WHERE eQTL = 1 AND %s
                    """ % (str(where)))
    results = cursor.fetchall()
    count_number = 0
    count_unique = 0
    set_snps = set()
    for res in results:
        print()
        count_number, set_snps, count_unique = close_eQTL(cursor, position_gene, int(res['pos_start']), int(res['pos_end']), int(res['chr']), where, count_number, set_snps, count_unique)
        # print(f"{res['pos_start']} - {res['pos_end']}")
    print('END')
    print(count_number)
    