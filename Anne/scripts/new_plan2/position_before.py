from Database import Database
import pandas as pd
# Python program to create
# sparse matrix using csr_matrix()
  
# Import required package
import numpy as np
from scipy.sparse import csr_matrix

def add_value(db):
    """

    """
    # add in_transcript
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_transcript` BOOLEAN DEFAULT(FALSE)
                    """)
    # add in_coding
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_coding` BOOLEAN DEFAULT(FALSE)
                    """)
    # add in_exon
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_exon` BOOLEAN DEFAULT(FALSE)
                    """)
    db.mydb_connection.commit()



def set_gene(db, index, row, chr):
    """

    """
    print(index)
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
    
    """
    # Make donor_list
    # List because double donors count
    donor_list = list()
    # Loop over list with snp IDs
    for snp_ID in snp_IDs:
        # Make donor_set
        # Set because duplicate snps from a donor don't count
        donor_set = set()
        # See how many snp_IDs there are in this table that are equal to the ID of the snp
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID
                        FROM 'donor_has_snp'
                        WHERE snp_ID = %s AND total_read_count > 0 AND mutant_allele_read_count > 0;
                        """ %
            (snp_ID)) #TODO total read count and mutation groter dan 0
        results = db.cursor.fetchall()
        for res in results:
            # Add donor_ID to the set
            donor_set.add(donor_dict[res['donor_ID']]) #res['donor_ID']
        # Extend the list with the donor_set
        donor_list.extend(donor_set)
    return donor_list


def close_to(db, gene, chr, start, end, position_out_gene, position_in_gene, before_after_gene, donor_dict, donor_list, gene_name_bef, gene_name_aft, sparseMatrix_before, sparseMatrix_after):
    # START
    print('START')
    gene_index_bef = gene_name_bef.index(gene)
    gene_name_bef[gene_index_bef] = f'{gene}:chr{chr}:{int(start)-position_out_gene}-{int(start)+position_in_gene}'
    db.cursor.execute("""
                    SELECT ID
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    """ %
        (str(chr), int(start)-position_out_gene, int(start)+position_in_gene))
    results_before = db.cursor.fetchall()
    # Make before_list
    before_list = list()
    print('append to list')
    for res_bef in results_before:
        # Add ID to before_list
        before_list.append(res_bef['ID'])
    # Call donor_list_before
    print('check donors')
    donor_list_before = check_snp_id(db, before_list, donor_dict)
    print('loop donors')
    for donor in donor_list_before:
        donor_index = donor_list.index(donor)
        sparseMatrix_before[donor_index, gene_index_bef] += 1

    # # END
    # db.cursor.execute("""
    #                 SELECT ID
    #                 FROM 'snp'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
    #                 """ %
    #     (str(chr), int(end)-position_in_gene, int(end)+position_out_gene))
    # results_after = db.cursor.fetchall()
    # # Make after_list
    # after_list = list()
    # for res_aft in results_after:
    #     # Add ID to after_list
    #     after_list.append(res_aft['ID'])
    # # Call donor_list_after
    # donor_list_after = check_snp_id(db, after_list, donor_dict)
    # # Write to file
    # before_after_gene.write(gene+'\t'+chr+'\t'+str(start)+'\t'+str(end)+'\t'+str(len(before_list))+'\t'+str(len(after_list))+'\t'+','.join(map(str,before_list))+'\t'+','.join(map(str,after_list))+'\t'+str(len(donor_list_before))+'\t'+str(len(donor_list_after))+'\t'+str(len(set(donor_list_before)))+'\t'+str(len(set(donor_list_after)))+'\t'+','.join(map(str,donor_list_before))+'\t'+','.join(map(str,donor_list_after))+'\n')

    return sparseMatrix_before, sparseMatrix_after, gene_name_bef, gene_name_aft

def write_sparse_matrix(sparseMatrix, gene_name_list, donor_list, save_path, pos):
    df = pd.DataFrame(data=sparseMatrix, columns=gene_name_list)
    df['donor_id'] = donor_list
    # df['cancer'] = donor_cancer_list
    df.set_index('donor_id')
    df.to_csv(f'{save_path}sparsematrix_{pos}.tsv.gz', sep="\t", index=True, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})

def loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, save_path, gene_name_bef, gene_name_aft, sparseMatrix_before, sparseMatrix_after):
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO
    # Make file
    before_after_gene = open(f'{save_path}gene_before_after_{position_out_gene}_{position_in_gene}.tsv', 'w')
    # Write header
    before_after_gene.write('gene\tchr\tstart_position\tend_position\tbefore\tafter\tbefore_list\tafter_list\tbefore_donor\tafter_donor\tbefore_donor_set\tafter_donor_set\tbefore_donor_list\tafter_donor_list\n')
    # Loop over genes in file
    for index, row in gene_df.iterrows():
        chr = row['chrom'].replace('chr', '')
        print(chr, row['txStart'], row['txEnd'])
        # Call set_gene
        set_gene(db, index, row, chr)
        # Call close_to
        sparseMatrix_before, sparseMatrix_after, gene_name_bef, gene_name_aft = close_to(db, row['#name'], chr, row['txStart'], row['txEnd'], position_out_gene, position_in_gene, before_after_gene, donor_dict, donor_list, gene_name_bef, gene_name_aft, sparseMatrix_before, sparseMatrix_after)
    # Close file
    before_after_gene.close()
    write_sparse_matrix(sparseMatrix_before, gene_name_bef, donor_list, save_path, 'bef')
    # write_sparse_matrix(sparseMatrix_after, gene_name_aft, donor_list, save_path, 'aft')
    

def get_donors(db):
    """
    """
    # DONOR
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    # Make dictionary with as ket ID and as value donor_ID
    donor_dict = dict()
    donor_set = set()
    for donor in donors:
        # Add to donor_dict
        donor_dict[donor['ID']] = donor['donor_ID']
        donor_set.add(donor['donor_ID'])
    return donor_dict, list(donor_set)

def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    # Path of the genes and there positions
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed" #snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
    # Save path
    save_path = "D:/Hanze_Groningen/STAGE/db/"
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    print(len(gene_df))
    gene_df['proteinID'].replace('', np.nan, inplace=True)
    gene_df.dropna(subset=['proteinID'], inplace=True)
    print(len(gene_df))
    gene_name_list = gene_df['#name'].tolist()
    # Positions before and after gene
    """
    Hoe Okosun et al. non-coding mutations in promotors aan gen/transcript aanwees: 
    * Regions of  -2000bp - 250bp (5' UTRs if applicable) from the transcription starting sites (TSS) 
    for each transcript were screened. For transcripts from the same gene that share the same promoter 
    mutation profiles, only one representative transcript was selected.
    """
    position_out_gene = 2000
    position_in_gene = 250
    # Call add_value
    # add_value(db)
    print('set GENE')
    # Call get_donors
    donor_dict, donor_list = get_donors(db)
    # Creating a len(donor_list) * len(gene_name_list) sparse matrix
    sparseMatrix_before = csr_matrix((len(donor_list), len(gene_name_list)), 
                            dtype = np.int8).toarray()
    sparseMatrix_after = csr_matrix((len(donor_list), len(gene_name_list)), 
                            dtype = np.int8).toarray()
    print('loop over genes')
    # Call loop_over_genes
    loop_over_genes(db, gene_df, position_out_gene, position_in_gene, donor_dict, donor_list, save_path, gene_name_list, gene_name_list.copy(), sparseMatrix_before, sparseMatrix_after)
    print('CLOSE')
    # Close database connection
    db.close()
    

if __name__ == '__main__':
    main()