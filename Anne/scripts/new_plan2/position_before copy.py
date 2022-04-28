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
    # Committing the current transactions
    db.mydb_connection.commit()


def set_gene(db, row, chr):
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
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(chr), int(row['hg19.knownGene.txStart']), int(row['hg19.knownGene.txEnd'])))
    # Update in_coding
    db.cursor.execute(
        """UPDATE snp
            SET in_coding = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
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
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(chr), int(exon_start[i]), int(exon_end[i])))
    # Committing the current transactions
    db.mydb_connection.commit()



def loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list,
                    save_path, gene_name_bef, gene_name_aft, sparse_matrix_before_region, sparse_matrix_after_region,
                    sparse_matrix_before_overall, sparse_matrix_after_overall, filter_num, with_type):
    """
    Loop over the gene data frame.
    :param db:                  The database object
    :param gene_df:             The data frame with all genes (columns: #name, chrom, strand, txStart, txEnd, cdsStart,
                                cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param position_out_gene:   Region before the start position of a gene or after the stop position of a gene
    :param position_in_gene:    Region after the start position of a gene or before the stop position of a gene
    :param donor_dict:          A dictionary with as key the automatically generated donor ID and as value the donor
                                IDs that are used in the research.
    :param donor_list:          List of donor names (to be used later as rows in the sparse matrix)
    :param donor_cancer_list:   List of cancers. This list has the same order as donor_list.
    :param save_path:           Path to save files
    :param gene_name_bef:       List of donor names before (to be used later as rows in the sparse matrix)
    :param gene_name_aft:       List of donor names after (to be used later as rows in the sparse matrix)
    :param sparse_matrix_before_region:A matrix which contains very few non-zero elements. It contains the counts of a
                                specific region (before gene)-donor combination.
    :param sparse_matrix_after_region: A matrix which contains very few non-zero elements. It contains the counts of a
                                specific region (after gene)-donor combination.
    :param sparse_matrix_before_overall:
    :param sparse_matrix_after_overall:
    :return:

    """
    # TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    # TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    # TODO
    # The header for the files before_gene_file and after_gene_file
    header_file = 'filter\tgene\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all' \
                  '\tdonor_count\tcancer_count\n'
    # Make file
    before_gene_file = open(f'{save_path}gene_before_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    before_gene_file.write(header_file)
    # Make file
    after_gene_file = open(f'{save_path}gene_after_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    after_gene_file.write(header_file)
    # TODO make nist
    total_read_before = [0] * len(donor_list)
    total_read_after = [0] * len(donor_list)
    # Loop over genes in file
    for index, row in gene_df.iterrows():
        # Remove 'chr' from the chromosome (chr1 --> 1)
        chr = row['hg19.knownGene.chrom'].replace('chr', '')
        print(chr, row['hg19.knownGene.txStart'], row['hg19.knownGene.txEnd'])
        # Call set_gene
        set_gene(db, row, chr)
        # Call close_to
        print('BEF')
        sparse_matrix_before_region, sparse_matrix_before_overall, gene_name_bef, total_read_before = close_to(db,
                                                                                                        row['hg19.kgXref.geneSymbol'],
                                                                                                        chr,
                                                                                                        row[
                                                                                                            'hg19.knownGene.txStart'] - position_out_gene,
                                                                                                        row[
                                                                                                            'hg19.knownGene.txStart'] + position_in_gene,
                                                                                                        before_gene_file,
                                                                                                        donor_dict,
                                                                                                        donor_list,
                                                                                                        gene_name_bef,
                                                                                                        sparse_matrix_before_region,
                                                                                                        sparse_matrix_before_overall,
                                                                                                        donor_cancer_list,
                                                                                                        total_read_before, filter_num, with_type)
        print('AFT')
        sparse_matrix_after_region, sparse_matrix_after_overall, gene_name_aft, total_read_after = close_to(db, row['hg19.kgXref.geneSymbol'],
                                                                                                      chr,
                                                                                                      row[
                                                                                                          'hg19.knownGene.txEnd'] - position_in_gene,
                                                                                                      row[
                                                                                                          'hg19.knownGene.txEnd'] + position_out_gene,
                                                                                                      after_gene_file,
                                                                                                      donor_dict,
                                                                                                      donor_list,
                                                                                                      gene_name_aft,
                                                                                                      sparse_matrix_after_region,
                                                                                                      sparse_matrix_after_overall,
                                                                                                      donor_cancer_list,
                                                                                                      total_read_after, filter_num, with_type)

    # Close file
    before_gene_file.close()
    after_gene_file.close()
    # Call write_sparse_matrix
    write_sparse_matrix(sparse_matrix_before_region, gene_name_bef, donor_list, save_path, 'bef_region',
                        donor_cancer_list, total_read_before)
    write_sparse_matrix(sparse_matrix_after_region, gene_name_aft, donor_list, save_path, 'aft_region',
                        donor_cancer_list, total_read_after)
    write_sparse_matrix(sparse_matrix_before_overall, gene_name_bef, donor_list, save_path, 'bef_overall',
                        donor_cancer_list, total_read_before)
    write_sparse_matrix(sparse_matrix_after_overall, gene_name_aft, donor_list, save_path, 'aft_overall',
                        donor_cancer_list, total_read_after)


def check_filter(db):
    db.cursor.execute("""
                    SELECT snp_ID, donor_ID, total_read_count, dosages, mutant_allele_read_count
                    FROM 'donor_has_snp'
                    WHERE total_read_count > 800 AND mutant_allele_read_count > 0;
                    """)
    results = db.cursor.fetchall()
    total_read_list = list()
    total_read_set = set()
    snp_set = set()
    for res in results:
        total_read_list.append(int(res['mutant_allele_read_count']))
        total_read_set.add(int(res['mutant_allele_read_count']))
        snp_set.add(int(res['snp_ID']))

    print(len(total_read_list))
    print(len(total_read_set))
    print(min(list(total_read_set)))
    print(max(list(total_read_set)))
    print('snp: ', len(snp_set))
    # count_read = dict(Counter(snp_set))
    # names = list(count_read.keys())
    # values = list(count_read.values())
    # plt.bar(range(len(count_read)), values, tick_label=names)
    # # plt.show()

    # # plt.hist(total_read_list, 50)#pd.Series(donor_list).hist()
    # plt.tight_layout()
    # pd.Series(total_read_list).plot.bar()
    # distribution = pd.Series(total_read_list).value_counts().sort_index()
    # print(distribution.head())
    
    # # distribution.plot.bar()
    # # distribution.head(200).plot.bar()
    # # distribution.iloc[:74].plot.bar()
    # # plt.tight_layout()
    # plt.savefig("D:/Hanze_Groningen/STAGE/DATAB/dist_var_0-74.png")
    # print(np.percentile(total_read_list, [25, 50, 75]))
    # print(pd.Series(total_read_list).describe())
    # print(pd.Series(total_read_list).quantile([0.25,0.5,0.75]))
    # print(pd.Series(total_read_list).quantile([0.05]))



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
    save_path = config['umap_path']
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
    filter_num = 33
    with_type = 'genes'
    # Call add_value
    # add_value(db)
    print('set GENE')
    # Call get_projects
    project_dict = db.get_projects()
    # Call get_donors
    donor_list, donor_dict, donor_cancer_list = db.get_donors(project_dict)
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    sparse_matrix_before_region = csr_matrix((len(donor_list), len(gene_name_list)),
                                             dtype=np.int8).toarray()
    sparse_matrix_after_region = csr_matrix((len(donor_list), len(gene_name_list)),
                                            dtype=np.int8).toarray()
    print('loop over genes')
    # Call loop_over_genes
    loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list,
                    save_path, gene_name_list, gene_name_list.copy(), sparse_matrix_before_region,
                    sparse_matrix_after_region, sparse_matrix_before_region.copy(), sparse_matrix_after_region.copy(), filter_num, with_type)
    print('CLOSE')
    # # Close database connection
    db.close()


if __name__ == '__main__':
    main()
