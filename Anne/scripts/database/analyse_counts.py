#!/usr/bin/env python3
import pandas as pd
import sys

import matplotlib.pyplot as plt
import numpy as np


from Database import Database
from germline_eQTL import place_germline
from gene import place_gene
from eQTL import place_eQTL
#TRUE = 1 and FALSE=0


def make_plot(xas, groups, title, in_trans, in_trans_germ, in_trans_somatic):
    #  names = ['False', 'True']
    # values = [results[0][1], results[1][1]]
    print(in_trans_germ[0][0])
    print(in_trans_germ[0][1])
    print(in_trans_germ[1][0])
    print(in_trans[0][0])
    print(in_trans[1][0])
    # intron = [in_trans_germ[0][0], in_trans_somatic[0][0]]
    # exon = [in_trans_germ[1][0], in_trans_somatic[1][0]]
    # r = [0, 1]
    
    # # Create brown bars
    # plt.bar(r, intron, color='#7f6d5f', edgecolor='white', width=1)
    # # Create green bars (middle), on top of the first ones
    # plt.bar(r, exon, bottom=intron, color='#557f2d', edgecolor='white', width=1)
    # # Custom X axis
    # plt.xticks(r, xas, fontweight='bold')
    # plt.title(title)
    # plt.xlabel('groups')
    # plt.ylabel('Counts')
    # plt.show()



def analyse_dbSNP_non_coding(db, cursor, gene_path, where, eQTL_path):
    # ##########
    # # dbSNP
    # ##########
    # print('dbSNP')
    # # Count germline vs somatic (dbsnp)
    # db.count_values('germline', 'snp', where, 'Germline (True) vs somatic (False) - dbSNP')
    # print('place_germline')
    # # Count mutations close to germline
    # place_germline(100, cursor, where)

    ##########
    # coding / non-coding
    ##########
    print('coding / non-coding')
    # Count in_transcript
    in_trans = db.count_values('in_transcript', 'snp', where, 'in transcript')
    # Count in_coding
    in_coding = db.count_values('in_coding', 'snp', where, 'in CDS (coding)')
    # Count in_exon
    in_exon = db.count_values('in_exon', 'snp', where, 'in exon')

    # ##########
    # # dbSNP & coding / non-coding
    # ##########
    # print('dbSNP & coding / non-coding')
    # value_germline = 1 #true
    # where_germline = f"seq_strategy = 'WGS' AND germline = {value_germline} "    
    # # Count in_transcript and germline == TRUE
    # in_trans_germ = db.count_values('in_transcript', 'snp', where_germline, 'in transcript germline')
    # # Count in_coding and germline == TRUE
    # # in_coding_germ = db.count_values('in_coding', 'snp', where_germline, 'in CDS (coding) germline')
    # # Count in_exon and germline == TRUE
    # # in_exon_germ = db.count_values('in_exon', 'snp', where_germline, 'in exon germline')
    # value_germline_FALSE = 0 
    # where_germline_FALSE = f"seq_strategy = 'WGS' AND germline = {value_germline_FALSE} "    
    # # Count in_transcript and germline == TRUE
    # in_trans_somatic = db.count_values('in_transcript', 'snp', where_germline_FALSE, 'in transcript germline')
    # # Count in_coding and germline == TRUE
    # # in_coding_somatic = db.count_values('in_coding', 'snp', where_germline_FALSE, 'in CDS (coding) germline')
    # # Count in_exon and germline == TRUE
    # # in_exon_somatic = db.count_values('in_exon', 'snp', where_germline_FALSE, 'in exon germline')
    # make_plot(['coding', 'non coding'], ['germline', 'somatic'], 'test', in_trans, in_trans_germ, in_trans_somatic)
    # # print('place_gene')
    # # # value moet 0 (false zijn), want het mag niet binnen een gen liggen
    # # valueGene = 0
    # # # Count mutations close to in_transcript
    # # place_gene(gene_path, 100, db.cursor, 'in_transcript', where_germline+f'AND in_transcript = {valueGene}', 'txStart', 'txEnd')
    # # place_gene(gene_path, 100, db.cursor, 'in_coding', where_germline+f'AND in_coding = {valueGene}', 'cdsStart', 'cdsEnd')
    # # place_gene(gene_path, 100, db.cursor, 'in_exon', where_germline+f'AND in_transcript = {valueGene}', 'exonStarts', 'exonEnds')

    # ##########
    # # eQTL
    # ##########
    # # mutation = eQTL
    # # db.count_values('eQTL', 'snp', where, 'eQTL')
    # # print('-------1--------')
    # # place_germline(100, db.cursor, where, 'eQTL')
    # # print('-------2--------')
    # # place_eQTL(eQTL_path, 100, cursor, where)

    


def main():
    # '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene_long.db'
    # db_path='/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/Database_internship_gene.db'
    db_path="D:/Hanze_Groningen/STAGE/DATAB/Database_internship_gene_long_NEW2.0 - kopie (3).db"
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed"
    eQTL_path = "C:/Users/Anne_/Downloads/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    db = Database(db_path) #sys.argv[1]
    
    where = "seq_strategy = 'WGS'"
    analyse_dbSNP_non_coding(db, db.cursor, gene_path, where, eQTL_path)
    db.close()
    
    


if __name__ == '__main__':
    main()
