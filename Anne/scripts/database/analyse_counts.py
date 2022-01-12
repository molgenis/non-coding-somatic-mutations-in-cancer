#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database
#TRUE = 1 and FALSE=0


def count_values_tets(cursor, column, table, where):
    """

    :param column:
    :param table:
    :return:
    """
    # self.cursor.execute(f"""
    #                     SELECT {column}, COUNT(*)
    #                     FROM {table}
    #                     GROUP BY {column};
    #                     """)
    # results = self.cursor.fetchall()
    # print(f'---{column}')
    # for res in results:
    #     print(f'{res[0]} - {res[1]}')
    print(f'FILTER: {where}')
    cursor.execute(f"""
                        SELECT {column}, COUNT(*)
                        FROM {table}
                        WHERE {where}
                        GROUP BY {column};
                        """)
    results = cursor.fetchall()
    print(results)
    for res in results:
        print(f'{res[0]} - {res[1]}')
        print(res)

    print('\n\n')


def close_gene(row, position_gene, cursor, column, table, where, start, end, count_number, set_snps, count_unique):
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

    count_number += len(start_set)
    # letters in a but not in set_snps
    unique_start = list(start_set - set(set_snps))
    count_unique += len(unique_start)
    set_snps.update(unique_start)
 

    # count SNPs
    # cursor.execute("""
    #                 SELECT count('%s')
    #                 FROM '%s'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND %s
    #                 """ %
    #     (str(column), str(table), str(row['chrom'].replace('chr', '')), 
    #     int(start)-position_gene, int(start), str(where)))
    # results = cursor.fetchall()
    # for res in results:
    #     print(res)
    #     count_number+=res[0]





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
    
    # cursor.execute("""
    #                 SELECT count('%s')
    #                 FROM '%s'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s AND %s
    #                 """ %
    #     (str(column), str(table), str(row['chrom'].replace('chr', '')), 
    #     int(end), int(end)+position_gene, str(where)))
    # results = cursor.fetchall()
    # for res in results:
    #     print(res)
    #     count_number+=res[0]

    return count_number, set_snps, count_unique



def place_gene(df_gene, position_gene, cursor, column, table, where, start, end):
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO

    gene_df = pd.read_csv(sys.argv[2], sep='\t')
    count_number = 0
    count_unique = 0
    set_snps = set()
    # CLOSE FROM GENES
    for index, row in gene_df.iterrows():
        if column == 'in_exon' and start == 'exonStarts' and end == 'exonEnds':
            exon_start = row[start].rstrip(',').split(',')
            exon_end = row[end].rstrip(',').split(',')
            print(f"COUNT - {row['exonCount']}")
            for i in range(int(row['exonCount'])):
                count_number = close_gene(row, position_gene, cursor, column, table, where, exon_start[i], exon_end[i], count_number, set_snps, count_unique)
        else:
            count_number = close_gene(row, position_gene, cursor, column, table, where, row[start], row[end], count_number, set_snps, count_unique)
            

    print(count_number)



def main():
    # '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long.db'
    # db_path='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene.db'
    db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene.db"
    db = Database(db_path) #sys.argv[1]
    
    COUNTER=0
    db.cursor.execute(f"""
                    SELECT *
                    FROM 'snp'
                    WHERE chr = '2';
                    """)
    results = db.cursor.fetchall()
    # print(results)
    for res in results:
        print(res['ID'])
        # COUNTER+=res[0]

    # cursor.execute(f"""
    #                 SELECT count(in_transcript)
    #                 FROM 'snp'
    #                 WHERE chr = '1' AND pos_start >= 1000000 AND pos_end <= 5000000 ;
    #                 """)
    # results = cursor.fetchall()
    # print(results)
    # for res in results:
    #     print(res)
    #     COUNTER+=res[0]

    # print(COUNTER)
    
    # where = "seq_strategy = 'WGS'"
    # # Count germline vs somatic (dbsnp)
    # # db.count_values('germline', 'snp')
    # # Count in_transcript
    # count_values_tets(cursor, 'in_transcript', 'snp', where)
    # # Count in_transcript
    # count_values_tets(cursor, 'in_coding', 'snp', where)
    # # Count in_transcript
    # count_values_tets(cursor, 'in_exon', 'snp', where)

    # where = "seq_strategy = 'WGS' AND "
    # # value=0 # false
    # value = 1 #true
    # # # Count in_transcript
    # # count_values_tets(cursor, 'in_transcript', 'snp', where+f'germline = {value}')
    # # count_values_tets(cursor, 'in_coding', 'snp', where+f'germline = {value}')
    # # count_values_tets(cursor, 'in_exon', 'snp', where+f'germline = {value}')
    # # value=1
    # # # Count in_transcript
    # # count_values_tets(cursor, 'in_transcript', 'snp', where+f'germline = {value}')
    # # count_values_tets(cursor, 'in_coding', 'snp', where+f'germline = {value}')
    # # count_values_tets(cursor, 'in_exon', 'snp', where+f'germline = {value}')


    # # value moet 0 (false zijn), want het mag niet binnen een gen liggen
    # place_gene('', 100, cursor, 'in_transcript', 'snp', where+f'in_transcript = {value}', 'txStart', 'txEnd')
    # place_gene('', 100, cursor, 'in_coding', 'snp', where+f'in_coding = {value}', 'cdsStart', 'cdsEnd')
    # place_gene('', 100, cursor, 'in_exon', 'snp', where+f'in_transcript = {value}', 'exonStarts', 'exonEnds')
    

    # # Count in_transcript
    # count_values_tets(cursor, 'in_transcript', 'snp', where+f'in_exon = {value}')
    # # Count in_transcript
    # count_values_tets(cursor, 'in_transcript', 'snp', where+f'in_exon = {value}')
    # # Count in_transcript
    # count_values_tets(cursor, 'in_transcript', 'snp', where+f'in_exon = {value}')

if __name__ == '__main__':
    main()
