import gzip
import pandas as pd
import os
import glob
# Python program to create
# sparse matrix using csr_matrix()
  
# Import required package
import numpy as np
import scipy.sparse
from scipy.sparse import csr_matrix


def set_donor_snp(project_file, set_donor, set_snp):
    header = project_file.readline().strip().split("\t")
    # Loop over the lines
    for num, line in enumerate(project_file):
        if (num %100000) == 0:
            print(num)
        # print(num)
        # Strip the line
        line = line.strip()
        # Split the line on tabs
        elems = line.split("\t")
        # Check if column 33 (=sequencing_strategy) is 'WGS'
        if elems[33] == 'WGS':
            donor_id = elems[1]
            set_donor.add(donor_id)
            chr = elems[8]
            pos = elems[9]
            ref = elems[15]
            alt = elems[16]
            snp_id = f'{chr}_{pos}_{ref}_{alt}'
            set_snp.add(snp_id)
    return set_donor, set_snp

def check_homo_hetero(total_read_count, mutant_allele_read_count):
    # Check if mutant_allele_read_count and total_read_count are empty
    # If these are not empty divide mutant_allele_read_count by total_read_count
    # and determine if they are homozygous or heterozygous.
    if mutant_allele_read_count != '' and total_read_count != '':
        if int(mutant_allele_read_count) != 0 or int(total_read_count) != 0:
            if int(mutant_allele_read_count) / int(total_read_count) >= 0.9:
                homo_hetero = 3 #'1/1'
            elif int(mutant_allele_read_count) / int(total_read_count) <= 0.1:
                homo_hetero = 1 #'0/0'
            else:
                homo_hetero = 2 #'0/1'
        else:
            homo_hetero = '.'

    # If mutant_allele_read_count and total_read_count are empty,
    # change both variables with a '.'.
    # Make homo_hetero a '.' for now. # TODO TODO
    else:
        homo_hetero = '.'
    return homo_hetero




def create_table(sparseMatrix, project_file, list_donor, list_snp):
    print('check')
    header = project_file.readline().strip().split("\t")
    snp_id_old = ''
    donor_id_old = ''
    
    # Loop over the lines
    for num, line in enumerate(project_file):
        if (num %100000) == 0:
            print(num)
        
        # print(num)
        # Strip the line
        line = line.strip()
        # Split the line on tabs
        elems = line.split("\t")
        # Check if column 33 (=sequencing_strategy) is 'WGS'
        if elems[33] == 'WGS':
            # icgc_mutation_id        icgc_donor_id   project_code 
            # Get the needed information out of de elems
            donor_id = elems[1]
            donor_index = list_donor.index(donor_id)
            # print(f'{donor_id} - {donor_index}')
            chr = elems[8]
            pos = elems[9]
            ref = elems[15]
            alt = elems[16]
            snp_id = f'{chr}_{pos}_{ref}_{alt}'
            if (snp_id_old != snp_id) and (donor_id_old != donor_id):
                snp_index = list_snp.index(snp_id)
                # print(f'{snp_id} - {snp_index}')

                total_read_count = elems[19]
                mutant_allele_read_count = elems[20]
                homo_hetero = check_homo_hetero(total_read_count, mutant_allele_read_count)
                if homo_hetero != '.':
                    # print(homo_hetero)
                    sparseMatrix[snp_index, donor_index] = homo_hetero
                # if project_code_old != project_code:
            snp_id_old = snp_id
            donor_id_old = donor_id
    return sparseMatrix

def create_vcf(sparseMatrix, list_donor, list_snp, name_vcf):
    # print('-----')
    df = pd.DataFrame(data=sparseMatrix, columns=list_donor)
    df['snp_id'] = list_snp
    df[['#CHROM', 'POS', 'REF', 'ALT']] = df['snp_id'].str.split('_', expand=True)
    df[['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT']] = '.'
    order_vcf = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + list_donor
    # print(order_vcf[:14])
    df = df[order_vcf]
    # print(df)
    df.to_csv(name_vcf, sep="\t", index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
            


def main():
    # The path to the data
    path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/'
    # The path where the new data should be stored
    out_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/vcf/'

    cancer_types_path = '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/self_made/site.csv' #"D:/Hanze_Groningen/STAGE/NEW PLAN/site.csv"
    cancer_types = pd.read_csv(cancer_types_path, sep=';')
    
    grouped = cancer_types.groupby(['cancer'])
    for cancer, projects in grouped:
        name_vcf = f'{out_path}{cancer}.vcf.gz'
        print(cancer)
        set_donor = set()
        set_snp = set()
        for project in projects['project_ID']:
            file = glob.glob(f"{path}download*.{project}.tsv.gz")
            # Open and unzip file
            project_file = gzip.open(file[0], 'rt')
            set_donor, set_snp = set_donor_snp(project_file, set_donor, set_snp)
        print(f'len donor: {len(set_donor)}')
        print(f'len snp: {len(set_snp)}')
        # Creating a len(list(set_snp)) * len(list(set_donor)) sparse matrix
        sparseMatrix = csr_matrix((len(list(set_snp)), len(list(set_donor))), 
                                dtype = np.int8).toarray()
        for project in projects['project_ID']:
            print(project)
            file = glob.glob(f"{path}download*.{project}.tsv.gz")
            # Open and unzip file
            project_file = gzip.open(file[0], 'rt')
            sparseMatrix = create_table(sparseMatrix, project_file, list(set_donor), list(set_snp))
        create_vcf(sparseMatrix, list(set_donor), list(set_snp), name_vcf)


if __name__ == '__main__':
    main()