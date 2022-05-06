from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix
# from matplotlib import pyplot as plt
from search_snps_between import close_to, write_sparse_matrix
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def add_value(db):
    """
    Adds values (in_transcript, in_coding, and in_exon) to the database (table snp).
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
    config = get_config()
    # Path of the database
    path_db = config['database'] #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" 
    # Database connection
    db = Database(path_db)
    
    # check_filter(db)
    # Path of the genes and there positions
    gene_path = config['all_genes'] # snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
    # Path to save files
    # save_path = config['umap_path']
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    # print(len(gene_df))
    # # Replace all empty values with NaN in the column proteinID
    # gene_df['proteinID'].replace('', np.nan, inplace=True)
    # # Drop all NaN values (in column proteinID)
    # gene_df.dropna(subset=['proteinID'], inplace=True)
    # print(len(gene_df))
    gene_name_list = gene_df['hg19.kgXref.geneSymbol'].tolist()
    print(len(gene_name_list))
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
    #https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-247
    #https://www.nature.com/articles/nature07517#Sec6
    # SNP discovery increases with increasing depth: essentially all homozygous positions are detected at 15×, whereas heterozygous positions accumulate more gradually to 33× (Fig. 5a). 
    # filter_num = 33
    # with_type = 'genes'
    # Call add_value
    add_value(db)
    print('set GENE')
    # Call get_projects
    # project_dict = db.get_projects()
    # Call get_donors
    # donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    print('loop over genes')
    # Call loop_over_genes
    loop_over_genes(db, gene_df, position_out_gene, position_in_gene)
    print('CLOSE')
    # # Close database connection
    db.close()


if __name__ == '__main__':
    main()
