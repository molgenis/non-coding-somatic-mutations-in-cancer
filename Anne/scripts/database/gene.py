#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database
#TRUE = 1 and FALSE=0


def close_gene(row, position_gene, cursor, where, start, end, count_number, set_snps, count_unique):
    # list(set(a) & set(b)) OR set(a).intersection(b) # overlap tussen twee lijsten
    # list(set(a() ^ set(b)) OR list(set(a).symmetric_difference(set(b)))
    """
    https://stackoverflow.com/questions/28444561/get-only-unique-elements-from-two-lists
    >>> a = set('abracadabra')
    >>> b = set('alacazam')
    >>> a                                  # unique letters in a
    {'a', 'r', 'b', 'c', 'd'}
    >>> a - b                              # letters in a but not in b
    {'r', 'd', 'b'}
    >>> a | b                              # letters in a or b or both
    {'a', 'c', 'r', 'd', 'b', 'm', 'z', 'l'}
    >>> a & b                              # letters in both a and b
    {'a', 'c'}
    >>> a ^ b                              # letters in a or b but not both
    {'r', 'd', 'b', 'm', 'z', 'l'}
    """
    

    # START
    cursor.execute("""
                    SELECT *
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND %s
                    """ %
        (str(row['chrom'].replace('chr', '')), int(start)-position_gene, int(start), str(where)))
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
        (str(row['chrom'].replace('chr', '')), int(end), int(end)+position_gene, str(where)))
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
def place_gene(gene_path, position_gene, cursor, column, where, start, end):
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO
    gene_df = pd.read_csv(gene_path, sep='\t')
    count_number = 0
    count_unique = 0
    set_snps = set()
    # CLOSE FROM GENES
    for index, row in gene_df.iterrows():
        if column == 'in_exon' and start == 'exonStarts' and end == 'exonEnds':
            exon_start = row[start].rstrip(',').split(',')
            exon_end = row[end].rstrip(',').split(',')
            for i in range(int(row['exonCount'])):
                count_number, set_snps, count_unique = close_gene(row, position_gene, cursor, where, exon_start[i], exon_end[i], count_number, set_snps, count_unique)
        else:
            count_number, set_snps, count_unique = close_gene(row, position_gene, cursor, where, row[start], row[end], count_number, set_snps, count_unique)
    print('END')
    print(count_number)