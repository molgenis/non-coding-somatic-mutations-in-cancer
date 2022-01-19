#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database
from germline import place_germline
from gene import place_gene
from analyse_eQTL import place_eQTL
#TRUE = 1 and FALSE=0

def analyse_dbSNP_non_coding(db, cursor, gene_path, where):
    # ##########
    # # dbSNP
    # ##########
    # print('dbSNP')
    # # Count germline vs somatic (dbsnp)
    # db.count_values('germline', 'snp', where)
    # print('place_germline')
    # # Count mutations close to germline
    # place_germline(100, cursor, where)

    # ##########
    # # coding / non-coding
    # ##########
    # print('coding / non-coding')
    # # Count in_transcript
    # db.count_values('in_transcript', 'snp', where)
    # # Count in_coding
    # db.count_values('in_coding', 'snp', where)
    # # Count in_exon
    # db.count_values('in_exon', 'snp', where)

    # ##########
    # # dbSNP & coding / non-coding
    # ##########
    # print('dbSNP & coding / non-coding')
    # value_germline = 1 #true
    # where_germline = f"seq_strategy = 'WGS' AND germline = {value_germline} "    
    # # Count in_transcript and germline == TRUE
    # db.count_values('in_transcript', 'snp', where_germline)
    # # Count in_coding and germline == TRUE
    # db.count_values('in_coding', 'snp', where_germline)
    # # Count in_exon and germline == TRUE
    # db.count_values('in_exon', 'snp', where_germline)

    # print('place_gene')
    # # value moet 0 (false zijn), want het mag niet binnen een gen liggen
    # valueGene = 0
    # # Count mutations close to in_transcript
    # place_gene(gene_path, 100, db.cursor, 'in_transcript', where_germline+f'AND in_transcript = {valueGene}', 'txStart', 'txEnd')
    # place_gene(gene_path, 100, db.cursor, 'in_coding', where_germline+f'AND in_coding = {valueGene}', 'cdsStart', 'cdsEnd')
    # place_gene(gene_path, 100, db.cursor, 'in_exon', where_germline+f'AND in_transcript = {valueGene}', 'exonStarts', 'exonEnds')

    ##########
    # eQTL
    ##########
    # mutation = eQTL
    db.count_values('eQTL', 'snp', where)
    place_germline(100, db.cursor, where, 'eQTL')
    


def main():
    # '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long.db'
    # db_path='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene.db'
    db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    db = Database(db_path) #sys.argv[1]
    
    where = "seq_strategy = 'WGS'"
    analyse_dbSNP_non_coding(db, db.cursor, gene_path, where)
    db.close()
    
    


if __name__ == '__main__':
    main()
