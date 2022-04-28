from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix
from matplotlib import pyplot as plt
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


def check_snp_id(db, snp_ids, donor_dict, donor_list, gene_index, sparse_matrix_overall, sparse_matrix_region,
                 total_read, filter_num):
    """
    Checks for each snp which donors have this snp. 
    For this snp, these donors must meet total_read_count > 0 and mutant_allele_read_count > 0.
    These donors are then added to a list. 
    This list will eventually contain the donors (a donor can be more frequent) who contain snps in a certain region.
    :param db:                The database object
    :param snp_ids:           A list of snp IDs that fall within a particular region.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :return: donor_list_snp : List of donors (research donor IDs) who contain snps within a particular region.
    """
    #
    donor_read_count = dict()
    # Make donor_list_snp. List because double donors count
    donor_list_snp = list()
    # Loop over list with snp IDs
    for snp_id in snp_ids:
        # Make donor_set. Set because duplicate snps from a donor don't count
        donor_set = set()
        # See how many snp_IDs there are in this table that are equal to the ID of the snp and total_read_count > 0
        # AND mutant_allele_read_count > 0
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID, total_read_count, dosages
                        FROM 'donor_has_snp'
                        WHERE snp_ID = %s AND total_read_count > %s AND mutant_allele_read_count > 0;
                        """ % (snp_id, filter_num))
        results = db.cursor.fetchall()
        for res in results:
            # Add donor_ID to the set
            donor_set.add(donor_dict[res['donor_ID']])
            if donor_dict[res['donor_ID']] in donor_read_count:
                donor_read_count[donor_dict[res['donor_ID']]][0] += res['total_read_count']
                donor_read_count[donor_dict[res['donor_ID']]][1] += res['dosages']
                donor_read_count[donor_dict[res['donor_ID']]][2] += 1
            else:
                donor_read_count[donor_dict[res['donor_ID']]] = [res['total_read_count'], res['dosages'], 1]
        # Extend the list with the donor_set
        donor_list_snp.extend(donor_set)
    # Loop over dict: donor_read_count
    for key, value in donor_read_count.items():
        donor_index = donor_list.index(key)
        sparse_matrix_overall[donor_index, gene_index] = value[2]
        # donor_read_count[key] = value[2] / (value[1] / value[0])
        sparse_matrix_region[donor_index, gene_index] = value[2] #/ (value[1] / value[0])
        if (value[2] / (value[1] / value[0])) <= 0:
            print(f'{value[2]} / ({value[1]} / {value[0]})')
        total_read[donor_index] += value[0]

    return donor_list_snp, sparse_matrix_overall, sparse_matrix_region, total_read

def filter_tot_read(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, sparse_matrix_region,
             sparse_matrix_overall, donor_cancer_list, total_read, snp_id_list, gene_index, filter_num):
    """

    :param db:                The database object
    :param gene:              The name of the gene
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param gene_file:         The file after which the genes with their specific region before or after
                              a gene are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param sparse_matrix_region: A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param sparse_matrix_overall:
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :param total_read:
    :param snp_id_list:
    :param gene_index:
    :param filter_num:           
    :return: donor_list_snp : List of donors (research donor IDs) who contain snps within a particular region.
    """
    if len(snp_id_list) > 0:
        # Call donor_list_before
        donor_list_snp, sparse_matrix_overall, sparse_matrix_region, total_read = check_snp_id(db, snp_id_list, donor_dict,
                                                                                            donor_list, gene_index,
                                                                                            sparse_matrix_overall,
                                                                                            sparse_matrix_region,
                                                                                            total_read, filter_num)
        # Make cancer_list
        cancer_list = list()
        for donor in donor_list_snp:
            donor_index = donor_list.index(donor)
            # Adds +1 to the sparse_matrix at the position of donor region
            # sparse_matrix_overall[donor_index, gene_index] += 1
            cancer_list.append(donor_cancer_list[donor_index])

        # Creates a dictionary from the list with as key the name (cancer type or donor ID) and
        # as value how often that name occurs in the list.
        cancer_count = dict(Counter(cancer_list))
        donor_count = dict(Counter(donor_list_snp))

        # Write to file
        gene_file.write(str(filter_num) + '\t'+gene + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_id_list)) + '\t' + ','.join(map(str, snp_id_list)) + '\t' + str(len(donor_list_snp)) + '\t' + str(
            donor_count) + '\t' + str(cancer_count) + '\n')
        # str(len(set(donor_list_snp)))+'\t'+','.join(map(str,donor_list_snp))+'\t'
    else:
        # Write to file
        gene_file.write(str(filter_num) + '\t'+gene + '\t' + chr + '\t' + str(start_pos) + '\t' + str(end_pos) + '\t' + str(
            len(snp_id_list)) + '\t-\t-\t-\t-\n')
    return sparse_matrix_region, sparse_matrix_overall, total_read



def close_to(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, gene_name_list, sparse_matrix_region,
             sparse_matrix_overall, donor_cancer_list, total_read):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (gene_file).
    :param db:                The database object
    :param gene:              The name of the gene
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param gene_file:         The file after which the genes with their specific region before or after
                              a gene are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor
                              IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param sparse_matrix_region:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return: sparse_matrix_region:   A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
             gene_name_list:  List of gene names (to be used later as columns in the sparse matrix)
    """
    # Get index of the gene name in the list
    gene_index = gene_name_list.index(gene)
    # Replace the name with gene_name:chromosoom:start_pos-end_pos
    gene_name_list[gene_index] = f'{gene}:chr{chr}:{int(start_pos)}-{int(end_pos)}'
    snp_id_list = db.get_snps(chr, start_pos, end_pos)
    sparse_matrix_region, sparse_matrix_overall, total_read = filter_tot_read(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, sparse_matrix_region,
             sparse_matrix_overall, donor_cancer_list, total_read, snp_id_list, gene_index, 0)
    # sparse_matrix_region, sparse_matrix_overall, total_read = filter_tot_read(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, sparse_matrix_region,
    #          sparse_matrix_overall, donor_cancer_list, total_read, snp_id_list, gene_index, 41)
    # sparse_matrix_region, sparse_matrix_overall, total_read = filter_tot_read(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, sparse_matrix_region,
    #          sparse_matrix_overall, donor_cancer_list, total_read, snp_id_list, gene_index, 50)
    

    return sparse_matrix_region, sparse_matrix_overall, gene_name_list, total_read


def write_sparse_matrix(sparse_matrix, gene_name_list, donor_list, save_path, pos, donor_cancer_list, total_read):
    """
    Writes the sparse_matrix to a compressed (.gz) .tsv file
    :param sparse_matrix:     A matrix which contains very few non-zero elements. It contains the counts of a specific
                              region-donor combination.
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param save_path:         Path to save files
    :param pos:               A string indicating whether it is before or after a gene (bef/aft)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return:

    """
    # Creates a dataframe from the sparse_matrix with column names gene_name_list
    df = pd.DataFrame(data=sparse_matrix, columns=gene_name_list)
    # Add donor IDs as column
    df['donor_id'] = donor_list
    # Add cancer types as column
    df['cancer'] = donor_cancer_list
    # Add total read count as column
    df['total_read'] = total_read
    # Makes the donor_ids column the index of the data frame
    df.set_index('donor_id', inplace=True)
    # Write the dataframe to a compressed .tsv file
    df.to_csv(f'{save_path}sparsematrix_{pos}.tsv.gz', sep="\t", index=True, encoding='utf-8',
              compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})


def loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list,
                    save_path, gene_name_bef, gene_name_aft, sparse_matrix_before_region, sparse_matrix_after_region,
                    sparse_matrix_before_overall, sparse_matrix_after_overall):
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
    total_read = [0] * len(donor_list)
    # Loop over genes in file
    for index, row in gene_df.iterrows():
        # Remove 'chr' from the chromosome (chr1 --> 1)
        chr = row['hg19.knownGene.chrom'].replace('chr', '')
        print(chr, row['hg19.knownGene.txStart'], row['hg19.knownGene.txEnd'])
        # Call set_gene
        set_gene(db, row, chr)
        # Call close_to
        print('BEF')
        sparse_matrix_before_region, sparse_matrix_before_overall, gene_name_bef, total_read = close_to(db,
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
                                                                                                        total_read)
        print('AFT')
        sparse_matrix_after_region, sparse_matrix_after_overall, gene_name_aft, total_read = close_to(db, row['hg19.kgXref.geneSymbol'],
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
                                                                                                      total_read)

    # Close file
    before_gene_file.close()
    after_gene_file.close()
    # Call write_sparse_matrix
    write_sparse_matrix(sparse_matrix_before_region, gene_name_bef, donor_list, save_path, 'bef_region',
                        donor_cancer_list, total_read)
    write_sparse_matrix(sparse_matrix_after_region, gene_name_aft, donor_list, save_path, 'aft_region',
                        donor_cancer_list, total_read)
    write_sparse_matrix(sparse_matrix_before_overall, gene_name_bef, donor_list, save_path, 'bef_overall',
                        donor_cancer_list, total_read)
    write_sparse_matrix(sparse_matrix_after_overall, gene_name_aft, donor_list, save_path, 'aft_overall',
                        donor_cancer_list, total_read)


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
    path_db = config['database']  #"D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db"
    # Database connection
    db = Database(path_db)
    
    # check_filter(db)
    # Path of the genes and there positions
    gene_path = config['all_genes'] #'D:/Hanze_Groningen/STAGE/db/all_genes_new.tsv' # snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
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
                    sparse_matrix_after_region, sparse_matrix_before_region.copy(), sparse_matrix_after_region.copy())
    print('CLOSE')
    # # Close database connection
    db.close()


if __name__ == '__main__':
    main()
