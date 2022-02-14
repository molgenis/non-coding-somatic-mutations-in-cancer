import gzip
import pandas as pd
from OverviewPerRow import OverviewPerRow

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
            if project_code_old != project_code:
                print(project_code)
                if project_code_old != '':
                    # Calls add_donors
                    overview = add_donors(overview)
                    # Calls make_vcf_file
                    make_vcf_file(overview, name_vcf, project_code_old)
                project_code_old = project_code
                # Make overview object
                overview = OverviewPerRow()
            # Call set_snp (function in the file: OVERVIEW_SNP)
            overview.set_snp(chr, pos, ref, alt, donor_id, total_read_count,
                            mutant_allele_read_count, specimen_id, project_code)

            
    # return overview

def check_dict(vcf_dict, all_donors):
    """

    """
    # Check whether all snps contain the same number of donors
    for key, value in vcf_dict.items():
        if len(value) != (len(all_donors)):
            print(f'--------------{key}')


def add_donors(overview):
    """

    """
    # print(overview.get_vcf_dict())
    # print(len(overview.get_all_donors()))
    # Get vcf_dict (key: SNP_ID, value: dict --> (key: vcf columns or donor_id,
    #                                              value: information out of columns or format structure of donor))
    vcf_dict = overview.get_vcf_dict()
    # Get all_donors (key: donor_id,
    #                 value: list of format structure <-- but this is for when the donors hasn't that snp)
    all_donors = overview.get_all_donors()

    # Loop over keys/values of vcf_dict
    for SNP_ID, vcf_info in vcf_dict.items():
        # Loop over keys/values all_donors
        for donor_ID, format_info in all_donors.items():
            # Check if donor_id exists in vcf_info (value of vcf_dict)
            # If it does not appear, add the donor with the empty format structure
            if donor_ID not in vcf_info:
                vcf_dict[SNP_ID][donor_ID] = format_info

    check_dict(vcf_dict, all_donors)
    return overview


def make_vcf_file(overview, name_vcf, project_code):
    """

    """
    # Get vcf_dict (key: SNP_ID, value: dict --> (key: vcf columns or donor_id,
    #                                              value: information out of columns or format structure of donor))
    vcf_dict = overview.get_vcf_dict()
    # Make dataframe with as row the values out vcf_dict
    df = pd.DataFrame.from_dict(list(vcf_dict.values()), orient='columns').sort_values(
        ["CHROM", "POS"]).reset_index().drop('index', axis=1)
    # Write file to csv and gzip that file
    df.to_csv(f'{name_vcf}vcf_{project_code}.tsv.gz', sep='\t', index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})






def main():
    # Path to the file
    path_file = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/allfiles.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/ALL-US.tsv.gz"
    # File name vcf file
    name_vcf = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/" #all_vcf.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/test_vcf.tsv.gz"
    # Open and unzip file
    project_file = gzip.open(path_file, 'rt')
    # Calls read_file
    read_file(project_file, name_vcf)
    


if __name__ == '__main__':
    main()