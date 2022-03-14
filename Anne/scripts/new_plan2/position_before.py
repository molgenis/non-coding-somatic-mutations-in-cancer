from Database import Database
import pandas as pd

def add_value(db):
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_transcript` BOOLEAN DEFAULT(FALSE)
                    """)
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_coding` BOOLEAN DEFAULT(FALSE)
                    """)
    # add dosages
    db.cursor.execute(f"""
                    ALTER TABLE snp
                    ADD `in_exon` BOOLEAN DEFAULT(FALSE)
                    """)
    db.mydb_connection.commit()



def set_gene(db, index, row):
    print(index)
    # Update in_transcript
    db.cursor.execute(
        """UPDATE snp
            SET in_transcript = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(row['chrom'].replace('chr', '')), int(row['txStart']), int(row['txEnd'])))
    # Update in_coding
    db.cursor.execute(
        """UPDATE snp
            SET in_coding = TRUE
            WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s;""" %
        (str(row['chrom'].replace('chr', '')), int(row['cdsStart']), int(row['cdsEnd'])))
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
            (str(row['chrom'].replace('chr', '')), int(exon_start[i]), int(exon_end[i])))
    
    db.mydb_connection.commit()
    

def check_snp_id(db, snp_IDs, donor_dict):
    donor_list = list()
    for snp_ID in snp_IDs:
        donor_set = set()
        db.cursor.execute("""
                        SELECT snp_ID, donor_ID
                        FROM 'donor_has_snp'
                        WHERE snp_ID = %s
                        """ %
            (snp_ID))
        results = db.cursor.fetchall()
        for res in results:
            donor_set.add(donor_dict[res['donor_ID']]) #res['donor_ID']
        donor_list.extend(donor_set)
    print(donor_list)
    return donor_list


def close_to(db, gene, chr, start, end, position_gene, before_after_gene, donor_dict):
    #START
    # db.cursor.execute("""
    #                 SELECT COUNT(ID)
    #                 FROM 'snp'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
    #                 """ %
    #     (str(chr.replace('chr', '')), int(start)-position_gene, int(start)))
    db.cursor.execute("""
                    SELECT ID
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    """ %
        (str(chr.replace('chr', '')), int(start)-position_gene, int(start)))
    results = db.cursor.fetchall()
    before_list = list()
    for res in results:
        before_list.append(res['ID'])
    donor_list_before = check_snp_id(db, before_list, donor_dict)
    
    db.cursor.execute("""
                    SELECT ID
                    FROM 'snp'
                    WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
                    """ %
        (str(chr.replace('chr', '')), int(end), int(end)+position_gene))
    results = db.cursor.fetchall()
    after_list = list()
    for res in results:
        after_list.append(res['ID'])
    donor_list_after = check_snp_id(db, after_list, donor_dict)

    

    before_after_gene.write(gene+'\t'+chr+'\t'+str(len(before_list))+'\t'+str(len(after_list))+'\t'+str(start)+'\t'+str(end)+'\t'+','.join(map(str,before_list))+'\t'+','.join(map(str,after_list))+'\t'+str(len(donor_list_before))+'\t'+str(len(donor_list_after))+'\t'+','.join(map(str,donor_list_before))+'\t'+','.join(map(str,donor_list_after))+'\n')
    
        
        


    
    # #START
    # db.cursor.execute("""
    #                 SELECT *
    #                 FROM 'snp'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
    #                 """ %
    #     (str(chr), int(start)-position_gene, int(start)))
    # results = db.cursor.fetchall()
    # start_set = set()
    # for res in results:
    #     start_set.add(res['ID'])

    # count_number += len(list(start_set))
    # # letters in a but not in set_snps
    # unique_start = list(start_set - set(set_snps))
    # count_unique += len(unique_start)
    # set_snps.update(unique_start)


    # #END
    # db.cursor.execute("""
    #                 SELECT *
    #                 FROM 'snp'
    #                 WHERE chr = '%s' AND pos_start >= %s AND pos_end <= %s
    #                 """ %
    #     (str(chr), int(end), int(end)+position_gene))
    # results = db.cursor.fetchall()
    # end_set = set()
    # for res in results:
    #     end_set.add(res['ID'])

    # count_number += len(end_set)
    # # letters in a but not in set_snps
    # unique_end = list(end_set - set(set_snps))
    # count_unique += len(unique_end)
    # set_snps.update(unique_end)

    # if count_number > 0:
    #     print(count_number)
    #     # print(set_snps)
    #     print(count_unique)
    #     print('-------')



def loop_over_DOT(db, gene_df, position_gene, donor_dict):
    #TODO Nu kunnen genen wel overlappen en dat uiteindelijk het gen wel binnen een ander gen ligt
    #TODO Het kan nu overlappen, er kan nu twee keer hetzelfde gen mee worden genomen
    #TODO
    before_after_gene = open(f'D:/Hanze_Groningen/STAGE/db/gene_before_after_{position_gene}.tsv', 'w')
    before_after_gene.write('gene\tchr\tbefore\tafter\tstart_position\tend_position\tbefore_list\tafter_list\tbefore_donor\tafter_donor\tbefore_donor_list\tafter_donor_list\n')
    for index, row in gene_df.iterrows():
        print(row['chrom'], row['cdsStart'], row['cdsEnd'])
        set_gene(db, index, row)
        #START, STOP
        close_to(db, row['#name'], row['chrom'], row['cdsStart'], row['cdsEnd'], position_gene, before_after_gene, donor_dict)
    before_after_gene.close()
        

    # # Add to database
    # db.mydb_connection.commit()

def get_donors(db):
    # DONOR
    db.cursor.execute("""SELECT *
                            FROM donor""")
    donors = db.cursor.fetchall()
    donor_dict = dict()
    for donor in donors:
        donor_dict[donor['ID']] = donor['donor_ID']
    return donor_dict

def main():
    # Path of the database
    path_db = "D:/Hanze_Groningen/STAGE/DATAB/copydatabase_C.db" #/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/new_db/copydatabase_C.db
    # Database connection
    db = Database(path_db)
    # Path of the genes and there positions
    gene_path = "D:/Hanze_Groningen/STAGE/db/snp132_ucsc_hg19_checkGene.bed" #snp132_ucsc_hg19_checkGene - kopie.bed , snp132_ucsc_hg19_checkGene.bed
    # Read gene file
    gene_df = pd.read_csv(gene_path, sep='\t')
    position_gene = 2000
    # print('add value')
    # add_value(db)
    print('set GENE')
    donor_dict = get_donors(db)
    loop_over_DOT(db, gene_df, position_gene, donor_dict)
    print('CLOSE')
    db.close()
    

if __name__ == '__main__':
    main()