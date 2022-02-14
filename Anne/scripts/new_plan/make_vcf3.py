import gzip
import pandas as pd
import numpy as np


def read_file(project_file, name_vcf):
    """

    """
    project_code_old = ''
    # mutation_id_old = ''
    # donor_id_old = ''
    # Headers of the files
    header = project_file.readline().strip().split("\t")
    
    # Loop over the lines
    for line in project_file:
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
            project_code = elems[2]
            specimen_id = elems[3]
            chr = elems[8]
            pos = elems[9]
            ref = elems[15]
            alt = elems[16]
            total_read_count = elems[19]
            mutant_allele_read_count = elems[20]
            SNP_ID = f'{chr}_{pos}_{ref}_{alt}'
            if project_code_old != project_code:
                project_code_old = project_code
                add_df = pd.DataFrame()
                add_df[SNP_ID, donor_id] = '0'
                print(add_df)
            else:
                add_df[SNP_ID, donor_id] = '0'

            # if project_code_old != project_code:
    #             print(project_code)
    #             if project_code_old != '':
    #                 # Calls add_donors
    #                 overview = add_donors(overview)
    #                 # Calls make_vcf_file
    #                 make_vcf_file(overview, name_vcf, project_code_old)
    #             project_code_old = project_code
    #             # Make overview object
    #             overview = OverviewPerRow()
    #             print(overview.dict_SNP_ID)
    #         # Call set_snp (function in the file: OVERVIEW_SNP)
    #         overview.set_snp(chr, pos, ref, alt, donor_id, total_read_count,
    #                         mutant_allele_read_count, specimen_id, project_code)
    #  # Calls add_donors
    # overview = add_donors(overview)
    # # Calls make_vcf_file
    # make_vcf_file(overview, name_vcf, project_code_old)
            
    # return overview





def main():
    # Path to the file
    path_file = "D:/Hanze_Groningen/STAGE/NEW PLAN/B.tsv.gz" #"/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/allfiles.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/ALL-US.tsv.gz"
    # File name vcf file
    name_vcf = "D:/Hanze_Groningen/STAGE/NEW PLAN/00" #"/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/" #all_vcf.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/test_vcf.tsv.gz"
    # Open and unzip file
    project_file = gzip.open(path_file, 'rt')

    # cmd = 'grep -c ".*" D:/Hanze_Groningen/STAGE/NEW PLAN/ALL_AML.tsv.gz' #f'zcat {path_file} | wc -l'
    # number_lines = os.system(cmd)
    # print(number_lines)
    # Calls read_file
    read_file(project_file, name_vcf)
    


if __name__ == '__main__':
    main()