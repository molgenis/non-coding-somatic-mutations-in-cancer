#!/usr/bin/env python3
import pandas as pd
import io
import gzip
import os
import sys



def read_vcf(path):
    """
    Read vcf files.
    Comes from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
        Author: Daichi Narushima (成島 大智)
        GitName: dceoy
        File: read_vcf.py
    :param path: Path to the file
    :return:
    """
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        lines = [l for l in f if not l.startswith('##')]
    if len(lines) != 0:
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                   'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        )#.rename(columns={'#CHROM': 'CHROM'})
    else:
        return []




def main():
    path_file = sys.argv[1] #"D:/Hanze_Groningen/STAGE/VCF/02eea485-625d-412b-b45a-458d87b1a866.vcf.gz"
    basename = os.path.basename(path_file).split('.')[0]
    print(basename)
    # cmd = f"cat {path_file} | grep '^##' | sed 's/=/\,/g' | sed 's/#//g' > {path}header.txt"
    # os.system(cmd)
    vcf_file = read_vcf(path_file)
    if len(vcf_file) != 0:
        vcf_file.rename(columns={'NORMAL': f'NORMAL_{basename}', 'TUMOR': f'TUMOR_{basename}'}, inplace=True)
        # print(vcf_file)
        vcf_file.to_csv(f'{sys.argv[2]}nohead.vcf', sep="\t", index=False, encoding='utf-8')
                    #compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
    else:
        print('NOOOO')
    


if __name__ == '__main__':
    main()