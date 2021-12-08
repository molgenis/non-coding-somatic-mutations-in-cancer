#!/usr/bin/env python3
import pandas as pd
import io
import sys
import os

# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    """

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

path = sys.argv[1]
print(path)
# Get the basename of the file
basename = os.path.basename(path) #.split('.')[0]
# Read vcf file
df = read_vcf(path)#(sys.argv[1].strip())
remove_dbSNP = df[~df['ID'].str.contains("rs")]
# write a dataframe to tsv file
remove_dbSNP.to_csv(f'{sys.argv[2]}noHeader_{basename}', sep="\t", index=False, header=None)

