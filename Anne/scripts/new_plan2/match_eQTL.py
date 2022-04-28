import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config




def compare_eQTL_files(df_eQTL, df_strong_eQTL, config):
    print('START')
    df = pd.DataFrame(columns=df_eQTL.columns)
    len_df = len(df)
    for index, row in df_strong_eQTL.iterrows():
        print(index)
        select_eQTL = df_eQTL[(df_eQTL['SNP']==row[0]) & (df_eQTL['Gene']==row[1])]
        df = df.append(select_eQTL, ignore_index=True)
        if len_df == len(df):
            print(len_df)
            print(row)
            print('----')
        len_df = len(df)
    df.to_csv(config['strong_eqtl_path'], sep='\t', encoding='utf-8', index=False)




def main():
    config = get_config()
    path_eQTL = config['all_sig_eqtl']
    df_eQTL = pd.read_csv(path_eQTL, sep='\t')
    print('read 1 done')
    path_strong_eQTL = config['strong_eqtl_noInfo']
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t', header=None)
    print('read 2 done')
    compare_eQTL_files(df_eQTL, df_strong_eQTL, config)



if __name__ == '__main__':
    main()