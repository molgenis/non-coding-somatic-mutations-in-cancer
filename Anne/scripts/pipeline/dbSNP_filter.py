#!/usr/bin/env python3
import sys
import os

# Also takes the folder 1 higher, so that I can do the import after
# sys.path.append("..")
sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf



def main():
    path = sys.argv[1]
    # Get the basename of the file
    basename = os.path.basename(path)
    # Read vcf file
    df = read_vcf(path)
    remove_dbSNP = df[~df['ID'].str.contains("rs")]
    # Write a dataframe to tsv file
    remove_dbSNP.to_csv(f'{sys.argv[2]}noHeader_{basename}', sep="\t", index=False, header=None)


if __name__ == "__main__":
    main()
