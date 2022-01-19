#!/usr/bin/env python3
import pandas as pd
import sys


from Database import Database
from germline import place_germline
from gene import place_gene
#TRUE = 1 and FALSE=0


def main():
    # '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long.db'
    # db_path='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene.db'
    db_path="D:/Hanze_Groningen/STAGE/TEST_DEL/Database_internship_gene_long_NEW2.0 - kopie.db"
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    db = Database(db_path) #sys.argv[1]
    cursor = db.cursor
    
    where = "seq_strategy = 'WGS'"
    
    ##########
    # dbSNP
    ##########
    # Count germline vs somatic (dbsnp)
    db.count_values('germline', 'snp')
    # Count mutations close to germline
    place_germline(100, cursor, where)

    ##########
    # coding / non-coding
    ##########
    # Count in_transcript
    db.count_values('in_transcript', 'snp', where)
    # Count in_coding
    db.count_values('in_coding', 'snp', where)
    # Count in_exon
    db.count_values('in_exon', 'snp', where)

    ##########
    # dbSNP & coding / non-coding
    ##########
    value_germline = 1 #true
    where = "seq_strategy = 'WGS' AND germline = {value_germline}"    
    # Count in_transcript and germline == TRUE
    db.count_values('in_transcript', 'snp', where)
    # Count in_coding and germline == TRUE
    db.count_values('in_coding', 'snp', where)
    # Count in_exon and germline == TRUE
    db.count_values('in_exon', 'snp', where)

    # value moet 0 (false zijn), want het mag niet binnen een gen liggen
    valueGene = 0
    # Count mutations close to in_transcript
    place_gene(gene_path, 100, db.cursor, 'in_transcript', where+f'AND in_transcript = {valueGene}', 'txStart', 'txEnd')
    place_gene(gene_path, 100, cursor, 'in_coding', 'snp', where+f'AND in_coding = {valueGene}', 'cdsStart', 'cdsEnd')
    place_gene(gene_path, 100, cursor, 'in_exon', 'snp', where+f'AND in_transcript = {valueGene}', 'exonStarts', 'exonEnds')


if __name__ == '__main__':
    main()
