import pandas as pd
import os
import io
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config




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
    with open(path, 'r', encoding='utf-8') as f: #gzip.open(path, 'rt', encoding='utf-8')
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
    config = get_config()
    path = "C:/Users/Anne_/Downloads/0c6ccd67-2fd2-45d9-a147-d1b6eb6a0b52.vcf"
    df = read_vcf(path)
    print(df.head())
    print(df.columns[9:])
    df.replace(r'0\|0', 1, regex=True, inplace=True)
    df.replace(r'0\|1', 2, regex=True, inplace=True)
    df.replace(r'1\|0', 2, regex=True, inplace=True)
    df.replace(r'1\|1', 3, regex=True, inplace=True)
    df.replace(r'./.', 0, regex=True, inplace=True)
    print(df.head())
    print(set(df['NORMAL']))
    print(set(df['TUMOR']))

    # path = config['vcf_path'] 
    # path_files = f"{path}*.vcf.gz"
    # # Loop over all files in path that ends with .tsv
    # for filename in glob.glob(path_files):
    #     # Get the base name of the specified path
    #     basename = os.path.basename(filename)    
    #     all_freq_path = f'{config['allel_freq']}af_{basename}'
    #     all_freq_vcf = f'{config['vcf_allel_freq']}vcf_af_{basename}'
    #     df = pd.read_csv(filename, sep='\t', compression='gzip')
    #     calculate_all_freq(df, all_freq_path, all_freq_vcf)


if __name__ == '__main__':
    main()