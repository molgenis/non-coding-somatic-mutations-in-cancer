#!/usr/bin/env python3

# Imports
import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from Database import Database
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def add_value(db):
    """
    Adds values (in_transcript, in_coding, in_exon, before_gene and after_gene) to the database (table snp).
    :param db:  The database object
    :return:
    """
    # Add in_transcript
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_transcript` BOOLEAN DEFAULT(FALSE)
                    """)
    # Add in_coding
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_coding` BOOLEAN DEFAULT(FALSE)
                    """)
    # Add in_exon
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_exon` BOOLEAN DEFAULT(FALSE)
                    """)
    # Add before_gene
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `before_gene` BOOLEAN DEFAULT(FALSE)
                    """)
    # Add after_gene
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `after_gene` BOOLEAN DEFAULT(FALSE)
                    """)
    # Committing the current transactions
    db.mydb_connection.commit()


def set_gene(db, row, chr, position_out_gene, position_in_gene):
    """
    Makes the values (in_transcript, in_coding, and in_exon) correct, by checking whether a snp is within certain
    start and stop.
    :param db:  The database object
    :param row: One row out of the gene file (columns: #name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
                exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param chr: Chromosome number or letter (without chr)
    :param position_out_gene: Region before the start position of a gene or after the stop position of a gene
    :param position_in_gene: Region after the start position of a gene or before the stop position of a gene
    :return:
    """
    # Update in_transcript
    db.cursor.execute(
        """UPDATE snp
            SET in_transcript = TRUE
            WHERE chr = '%s' AND pos_start BETWEEN %s and %s;""" %
        (str(chr), int(row['hg19.knownGene.txStart']), int(row['hg19.knownGene.txEnd'])))
    # Update before_gene
    db.cursor.execute(
        """UPDATE snp
            SET before_gene = TRUE
            WHERE chr = '%s' AND pos_start BETWEEN %s and %s;""" %
        (str(chr), int(row['hg19.knownGene.txStart']-position_out_gene), int(row['hg19.knownGene.txStart']+position_in_gene)))
    # Update after_gene
    db.cursor.execute(
        """UPDATE snp
            SET after_gene = TRUE
            WHERE chr = '%s' AND pos_start BETWEEN %s and %s;""" %
        (str(chr), int(row['hg19.knownGene.txEnd'])-position_in_gene, int(row['hg19.knownGene.txEnd'])+position_out_gene))
    # Update in_coding
    db.cursor.execute(
        """UPDATE snp
            SET in_coding = TRUE
            WHERE chr = '%s' AND pos_start BETWEEN %s and %s;""" %
        (str(chr), int(row['hg19.knownGene.cdsStart']), int(row['hg19.knownGene.cdsEnd'])))
    # Get start and end of the exons
    exon_start = row['hg19.knownGene.exonStarts'].rstrip(',').split(',')
    exon_end = row['hg19.knownGene.exonEnds'].rstrip(',').split(',')
    # Loop over the exons start-end
    for i in range(int(row['hg19.knownGene.exonCount'])):
        # Update in_exon
        db.cursor.execute(
            """UPDATE snp
                SET in_exon = TRUE
                WHERE chr = '%s' AND pos_start BETWEEN %s and %s;""" %
            (str(chr), int(exon_start[i]), int(exon_end[i])))
    # Committing the current transactions
    db.mydb_connection.commit()



def loop_over_genes(db, gene_df, position_out_gene, position_in_gene):
    """
    Loop over the gene data frame.
    :param db:                  The database object
    :param gene_df:             The data frame with all genes (columns: #name, chrom, strand, txStart, txEnd, cdsStart,
                                cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param position_out_gene:   Region before the start position of a gene or after the stop position of a gene
    :param position_in_gene:    Region after the start position of a gene or before the stop position of a gene
    :return:

    """
   
    for index, row in gene_df.iterrows():
        # Remove 'chr' from the chromosome (chr1 --> 1)
        chr = row['hg19.knownGene.chrom'].replace('chr', '')
        print(chr, row['hg19.knownGene.txStart'], row['hg19.knownGene.txEnd'])
        # Call set_gene
        set_gene(db, row, chr, position_out_gene, position_in_gene)


def main():
    # Call get_config
    config = get_config('gearshift')
    # Path of the database
    path_db = config['database'] 
    # Database connection
    db = Database(path_db)
    # Path of the genes and there positions
    gene_path = config['all_genes']
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    """
    Okosun et al. 
    * Regions of  -2000bp - 250bp (5' UTRs if applicable) from the transcription starting sites (TSS) 
    for each transcript were screened. For transcripts from the same gene that share the same promoter 
    mutation profiles, only one representative transcript was selected.
    """
    # Region before the start position of a gene or after the stop position of a gene
    position_out_gene = 2000
    # Region after the start position of a gene or before the stop position of a gene
    position_in_gene = 250
    # Call add_value
    add_value(db)
    # Call loop_over_genes
    loop_over_genes(db, gene_df, position_out_gene, position_in_gene)
    # Close database connection
    db.close()


if __name__ == '__main__':
    main()
