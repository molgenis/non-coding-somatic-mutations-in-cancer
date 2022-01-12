#!/usr/bin/env python3
import pandas as pd
import io
import sys
import os

# Also takes the folder 1 higher, so that I can do the import after
sys.path.append("..")
from vcf_reader import read_vcf

path = sys.argv[1]
print(path)
# Get the basename of the file
basename = os.path.basename(path) #.split('.')[0]
# Read vcf file
df = read_vcf(path)#(sys.argv[1].strip())
remove_dbSNP = df[~df['ID'].str.contains("rs")]
# write a dataframe to tsv file
remove_dbSNP.to_csv(f'{sys.argv[2]}noHeader_{basename}', sep="\t", index=False, header=None)

