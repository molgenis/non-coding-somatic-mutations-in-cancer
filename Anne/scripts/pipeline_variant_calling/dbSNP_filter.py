#!/usr/bin/env python3

# Imports
import sys
import os
sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf

def main():
    """
    The SNPs that match the dbSNP are removed, because these might be already
    present germline mutations.
    """
    path = sys.argv[1]
    # Get the basename of the file
    basename = os.path.basename(path)
    # Read vcf file
    df = read_vcf(path)
    # The SNPs that match the dbSNP are removed, because these might be already present germline mutations.
    remove_dbSNP = df[~df['ID'].str.contains("rs")]
    # Write a dataframe to tsv file
    remove_dbSNP.to_csv(f'{sys.argv[2]}noHeader_{basename}', sep="\t", index=False, header=None)


if __name__ == "__main__":
    main()
