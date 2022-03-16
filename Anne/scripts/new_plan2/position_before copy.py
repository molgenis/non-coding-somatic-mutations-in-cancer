from Database import Database
import pandas as pd
from collections import Counter
# Python program to create
# sparse matrix using csr_matrix()
# Import required package
import numpy as np
from scipy.sparse import csr_matrix

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
    Makes the values (in_transcript, in_coding, and in_exon) correct, by checking whether a snp is within certain start and stop.
    :param db:  The database object
    :param row: One row out of the gene file (columns: #name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param chr: Chromosome number or letter (without chr)
    :return:
    """
    # Update in_transcript
    db.cursor.execute(
        """UPDATE snp
            SET in_transcript = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(chr), int(row['txStart']), int(row['txEnd'])))
    # Update in_coding
    db.cursor.execute(
        """UPDATE snp
            SET in_coding = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(chr), int(row['cdsStart']), int(row['cdsEnd'])))
    # Get start and end of the exons
    exon_start = row['exonStarts'].rstrip(',').split(',')
    exon_end = row['exonEnds'].rstrip(',').split(',')
    # Loop over the exons start-end
    for i in range(int(row['exonCount'])):
        # Update in_exon
        db.cursor.execute(
            """UPDATE snp
                SET in_exon = TRUE
                WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
            (str(chr), int(exon_start[i]), int(exon_end[i])))
    # Committing the current transactions
    db.mydb_connection.commit()
    

def check_snp_id(db, snp_IDs, donor_dict):
    """
    Checks for each snp which donors have this snp. 
    For this snp, these donors must meet total_read_count > 0 and mutant_allele_read_count > 0.
    These donors are then added to a list. 
    This list will eventually contain the donors (a donor can be more frequent) who contain snps in a certain region.
    :param db:                The database object
    :param snp_IDs:           A list of snp IDs that fall within a particular region.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor IDs that are used in the research.
    :return: donor_list_snp : List of donors (research donor IDs) who contain snps within a particular region.
    """
    # Make donor_list_snp. List because double donors count
    donor_list_snp = list()
    # Loop over list with snp IDs
    for snp_ID in snp_IDs:
        # Make donor_set. Set because duplicate snps from a donor don't count
        donor_set = set()
        # See how many snp_IDs there are in this table that are equal to the ID of the snp and total_read_count > 0 AND mutant_allele_read_count > 
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID
                        FROM 'donor_has_snp'
                        WHERE snp_ID = %s AND total_read_count > 0 AND mutant_allele_read_count > 0;
                        """ %
            (snp_ID))
        results = db.cursor.fetchall()
        for res in results:
            # Add donor_ID to the set
            donor_set.add(donor_dict[res['donor_ID']]) #res['donor_ID']
        # Extend the list with the donor_set
        donor_list_snp.extend(donor_set)
    return donor_list_snp


def close_to(db, gene, chr, start_pos, end_pos, gene_file, donor_dict, donor_list, gene_name_list, sparseMatrix, donor_cancer_list):
    """
    Selects all snps that occur in a specific region on a specific chromosome.
    And writes a line with certain information in the file (gene_file).
    :param db:                The database object
    :param gene:              The name of the gene
    :param chr:               Chromosome number or letter (without chr)
    :param start_pos:         The starting position or a particular region
    :param end_pos:           The stop position of a particular region
    :param gene_file:         The file after which the genes with their specific region before or after a gene are written.
    :param donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor IDs that are used in the research.
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param sparseMatrix:      A matrix which contains very few non-zero elements. It contains the counts of a specific region-donor combination.
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return: sparseMatrix:    A matrix which contains very few non-zero elements. It contains the counts of a specific region-donor combination.
             gene_name_list:  List of gene names (to be used later as columns in the sparse matrix)
    """
    # Get index of the gene name in the list
    gene_index = gene_name_list.index(gene)
    # Replace the name with gene_name:chromosoom:start_pos-end_pos
    gene_name_list[gene_index] = f'{gene}:chr{chr}:{int(start_pos)}-{int(end_pos)}'
    # Find all snps that are on a certain chromosome in a certain region (between a certain start and stop)
    db.cursor.execute("""
                    SELECT ID
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    """ %
        (str(chr), int(start_pos), int(end_pos)))
    results_before = db.cursor.fetchall()
    # Make snp_id_list
    snp_id_list = list()
    for res_bef in results_before:
        # Add ID to snp_id_list
        snp_id_list.append(res_bef['ID'])
    # Call donor_list_before
    donor_list_snp = check_snp_id(db, snp_id_list, donor_dict)
    # Make cancer_list
    cancer_list = list()
    for donor in donor_list_snp:
        donor_index = donor_list.index(donor)
        # Adds +1 to the sparseMatrix at the position of donor region
        sparseMatrix[donor_index, gene_index] += 1
        cancer_list.append(donor_cancer_list[donor_index])

    # Creates a dictionary from the list with as key the name (cancer type or donor ID) and as value how often that name occurs in the list.
    cancer_count = dict(Counter(cancer_list))
    donor_count = dict(Counter(donor_list_snp))

    # Write to file
    gene_file.write(gene+'\t'+chr+'\t'+str(start_pos)+'\t'+str(end_pos)+'\t'+str(len(snp_id_list))+'\t'+','.join(map(str,snp_id_list))+'\t'+str(len(donor_list_snp))+'\t'+str(donor_count)+'\t'+str(cancer_count)+'\n') #str(len(set(donor_list_snp)))+'\t'+','.join(map(str,donor_list_snp))+'\t'

    return sparseMatrix, gene_name_list


def write_sparse_matrix(sparseMatrix, gene_name_list, donor_list, save_path, pos, donor_cancer_list):
    """
    Writes the sparseMatrix to a compressed (.gz) .tsv file
    :param sparseMatrix:      A matrix which contains very few non-zero elements. It contains the counts of a specific region-donor combination.
    :param gene_name_list:    List of gene names (to be used later as columns in the sparse matrix)
    :param donor_list:        List of donor names (to be used later as rows in the sparse matrix)
    :param save_path:         Path to save files
    :param pos:               A string indicating whether it is before or after a gene (bef/aft)
    :param donor_cancer_list: List of cancers. This list has the same order as donor_list.
    :return:

    """
    # Creates a dataframe from the sparseMatrix with column names gene_name_list
    df = pd.DataFrame(data=sparseMatrix, columns=gene_name_list)
    # Add donor IDs as column
    df['donor_id'] = donor_list
    # Add cancer types as column
    df['cancer'] = donor_cancer_list
    # Makes the donor_ids column the index of the data frame
    df.set_index('donor_id', inplace = True)
    # Write the dataframe to a compressed .tsv file
    df.to_csv(f'{save_path}sparsematrix_{pos}.tsv.gz', sep="\t", index=True, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})



def loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list, save_path, gene_name_bef, gene_name_aft, sparseMatrix_before, sparseMatrix_after):
    """
    Loop over the gene data frame.
    :param db:                  The database object
    :param gene_df:             The data frame with all genes (columns: #name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, proteinID, alignID)
    :param position_out_gene:   Region before the start position of a gene or after the stop position of a gene
    :param position_in_gene:    Region after the start position of a gene or before the stop position of a gene
    :param donor_dict:          A dictionary with as key the automatically generated donor ID and as value the donor IDs that are used in the research.
    :param donor_list:          List of donor names (to be used later as rows in the sparse matrix)
    :param donor_cancer_list:   List of cancers. This list has the same order as donor_list.
    :param save_path:           Path to save files
    :param gene_name_bef:       List of donor names before (to be used later as rows in the sparse matrix)
    :param gene_name_aft:       List of donor names after (to be used later as rows in the sparse matrix)
    :param sparseMatrix_before: A matrix which contains very few non-zero elements. It contains the counts of a specific region (before gene)-donor combination.
    :param sparseMatrix_after:  A matrix which contains very few non-zero elements. It contains the counts of a specific region (after gene)-donor combination.
    :return:
    
    """
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO
    # The header for the files before_gene_file and after_gene_file
    header_file = 'gene\tchr\tstart_position_regio\tend_position_regio\t#snp_unique\tsnp_list\t#donors_all\tdonor_count\tcancer_count\n'
    # Make file
    before_gene_file = open(f'{save_path}gene_before_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    before_gene_file.write(header_file)
    # Make file
    after_gene_file = open(f'{save_path}gene_after_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    after_gene_file.write(header_file)
    # Loop over genes in file
    for index, row in gene_df.iterrows():
        # Remove 'chr' from the chromosome (chr1 --> 1)
        chr = row['chrom'].replace('chr', '')
        print(chr, row['txStart'], row['txEnd'])
        # Call set_gene
        set_gene(db, row, chr)
        # Call close_to
        print('BEF')
        sparseMatrix_before, gene_name_bef = close_to(db, row['#name'], chr, row['txStart']-position_out_gene, row['txStart']+position_in_gene, before_gene_file, donor_dict, donor_list, gene_name_bef, sparseMatrix_before, donor_cancer_list)
        print('AFT')
        sparseMatrix_after, gene_name_aft = close_to(db, row['#name'], chr, row['txEnd']-position_in_gene, row['txEnd']+position_out_gene, after_gene_file, donor_dict, donor_list, gene_name_aft, sparseMatrix_after, donor_cancer_list)

    # Close file
    before_gene_file.close()
    after_gene_file.close()
    # Call write_sparse_matrix
    write_sparse_matrix(sparseMatrix_before, gene_name_bef, donor_list, save_path, 'bef', donor_cancer_list)
    write_sparse_matrix(sparseMatrix_after, gene_name_aft, donor_list, save_path, 'aft', donor_cancer_list)
    

def get_projects(db):
    """
    Searches for all projects in the database and converts them into a specific dictionary. (key: ID, value: cancer)
    :param db:             The database object
    :return: project_dict: Dictionary with as key project ID (automatically generated) and as value the type of cancer that belongs to it.
    """
    # Select all projects in the database
    db.cursor.execute("""SELECT *
                            FROM project""")
    projects = db.cursor.fetchall()
    project_dict = dict()
    for proj in projects:
        # Add to project_dict as key the auto-generated ID and as value the type of cancer that belongs to that project
        project_dict[proj['ID']] = proj['cancer']
    return project_dict



def get_donors(db, project_dict):
    """
    Search for all donors in the database and make different things with them (lists and dictionaries).
    :param db:                  The database object
    :param project_dict:        Dictionary with as key project ID (automatically generated) and as value the type of cancer that belongs to it.
    :return: donor_list:        List of donor names (to be used later as rows in the sparse matrix)
             donor_dict:        A dictionary with as key the automatically generated donor ID and as value the donor IDs that are used in the research.
             donor_cancer_list: List of cancers. This list has the same order as donor_list.
    """
    # Select all donors in the database
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    donor_dict = dict()
    donor_list = list()
    donor_cancer_list = list()
    for donor in donors:
        # Add to donor_dict as key the auto-generated ID and as value the donor ID as used in the study
        donor_dict[donor['ID']] = donor['donor_ID']
        donor_list.append(donor['donor_ID'])
        donor_cancer_list.append(project_dict[donor['project_ID']])
    return donor_list, donor_dict, donor_cancer_list






def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    # Path of the genes and there positions
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed" #snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
    # Path to save files
    save_path = "D:/Hanze_Groningen/STAGE/db/"
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    print(len(gene_df))
    # Replace all empty values with NaN in the column proteinID
    gene_df['proteinID'].replace('', np.nan, inplace=True)
    # Drop all NaN values (in column proteinID)
    gene_df.dropna(subset=['proteinID'], inplace=True)
    print(len(gene_df))
    gene_name_list = gene_df['#name'].tolist()
    print(len(gene_name_list))
    """
    Hoe Okosun et al. non-coding mutations in promotors aan gen/transcript aanwees: 
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
    project_dict = get_projects(db)
    # Call get_donores
    donor_list, donor_dict, donor_cancer_list = get_donors(db, project_dict)    
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    sparseMatrix_before = csr_matrix((len(donor_list), len(gene_name_list)), 
                            dtype = np.int8).toarray()
    sparseMatrix_after = csr_matrix((len(donor_list), len(gene_name_list)), 
                            dtype = np.int8).toarray()
    print('loop over genes')
    # Call loop_over_genes
    loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, donor_cancer_list, save_path, gene_name_list, gene_name_list.copy(), sparseMatrix_before, sparseMatrix_after)
    print('CLOSE')
    # Close database connection
    db.close()
    

if __name__ == '__main__':
    main()