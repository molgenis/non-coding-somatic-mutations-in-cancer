#!/usr/bin/env python3

# Imports
import pandas as pd
import os
import glob
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def calculate_all_freq(df, all_freq_path, all_freq_vcf):
    """
    Calculate the AF and writes it to a vcf
    :param df:  
    :param all_freq_path:
    :param all_freq_vcf:
    :return:
    """
    format_index = list(df.columns).index('FORMAT') + 1
    all_freq_df = df.iloc[:, :format_index] 
    select_df = df.iloc[:, format_index:] 
    # Count how many 0,1,2 and NA there are
    all_freq_df['NA'] = select_df.eq(0).sum(axis=1)
    all_freq_df['0'] = select_df.eq(1).sum(axis=1)
    all_freq_df['1'] = select_df.eq(2).sum(axis=1)
    all_freq_df['2'] = select_df.eq(3).sum(axis=1)
    
    all_freq_df['sum_alt'] = (all_freq_df['0'] * 0) + (all_freq_df['1'] * 1) + (all_freq_df['2'] * 2)
    all_freq_df['donors'] = all_freq_df.loc[:,['0', '1', '2']].sum(axis = 1)
    # Calculate AF
    all_freq_df['alle_freq_alt'] = all_freq_df['sum_alt'] / (all_freq_df['donors'] * 2)
    all_freq_df['alle_freq_ref'] = 1 - all_freq_df['alle_freq_alt']
    all_freq_df['donors_nan'] = all_freq_df.loc[:,['NA', '0', '1', '2']].sum(axis = 1)
    all_freq_df['alle_freq_nan_alt'] = all_freq_df['sum_alt'] / (all_freq_df['donors_nan'] * 2)
    all_freq_df['alle_freq_nan_ref'] = 1 - all_freq_df['alle_freq_nan_alt']
    all_freq_df.to_csv(all_freq_path, sep="\t", index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
    df['INFO'] = 'AF=' + all_freq_df['alle_freq_alt'].astype(str) + ':AFN=' + all_freq_df['alle_freq_nan_alt'].astype(str) 
    df.to_csv(all_freq_vcf, sep="\t", index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})

    

def main():    
    # Call get_config
    config = get_config('gearshift')
    path = config['vcf_path'] 
    path_files = f"{path}*.vcf.gz"
    # Loop over all files in path that ends with .tsv
    for filename in glob.glob(path_files):
        # Get the base name of the specified path
        basename = os.path.basename(filename)    
        all_freq_path = f'{config["allel_freq"]}af_{basename}'
        all_freq_vcf = f'{config["vcf_allel_freq"]}vcf_af_{basename}'
        df = pd.read_csv(filename, sep='\t', compression='gzip')
        # Call calculate_all_freq
        calculate_all_freq(df, all_freq_path, all_freq_vcf)


if __name__ == '__main__':
    main()




