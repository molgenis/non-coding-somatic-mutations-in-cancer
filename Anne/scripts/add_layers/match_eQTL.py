#!/usr/bin/env python3
import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

def compare_eQTL_files(df_eQTL, df_strong_eQTL, config):
    """
    Appends the extra information from the file with all eQTLs to the file with only the strictest eQTLs per gene.
    :param df_eQTL: all significant eQTLs
    :param df_strong_eQTL: all strongest eQTLs per gene (no information further)
    :param config: dictionary of all paths
    :return:
    """
    # Make empty dataframe with the same column names as df_eQTL
    df = pd.DataFrame(columns=df_eQTL.columns)
    # Loop over df_strong_eQTL
    for index, row in df_strong_eQTL.iterrows():
        # Get the info of the strong eQTL
        select_eQTL = df_eQTL[(df_eQTL['SNP']==row[0]) & (df_eQTL['Gene']==row[1])]
        # Append to df
        df = df.append(select_eQTL, ignore_index=True)
    # Write df to file
    df.to_csv(config['strong_eqtl_path'], sep='\t', encoding='utf-8', index=False)


def main():
    config = get_config('gearshift')
    # Path to all significant eQTLs
    path_eQTL = config['all_sig_eqtl']
    # Read file
    df_eQTL = pd.read_csv(path_eQTL, sep='\t')
    # Path to all strongest eQTLs per gene (no information further)
    path_strong_eQTL = config['strong_eqtl_noInfo']
    # Read file
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t', header=None)
    # Call compare_eQTL_files
    compare_eQTL_files(df_eQTL, df_strong_eQTL, config)



if __name__ == '__main__':
    main()