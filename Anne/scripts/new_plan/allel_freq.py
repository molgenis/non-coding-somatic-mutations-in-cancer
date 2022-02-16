import gzip
import pandas as pd
import os
import glob



def calculate_all_freq(df, all_freq_path, nan_value):
    format_index = list(df.columns).index('FORMAT') + 1
    all_freq_df = df.iloc[:, :format_index] 
    select_df = df.iloc[:, format_index:] 
    if nan_value:    
        all_freq_df['NA'] = select_df.eq(0).sum(axis=1)
        all_freq_df['0'] = select_df.eq(1).sum(axis=1)
    else: 
        all_freq_df['0'] = select_df.eq(1).sum(axis=1) + select_df.eq(0).sum(axis=1)
    all_freq_df['1'] = select_df.eq(2).sum(axis=1)
    all_freq_df['2'] = select_df.eq(3).sum(axis=1)
    all_freq_df['sum'] = (all_freq_df['0'] * 0) + (all_freq_df['1'] * 1) + (all_freq_df['2'] * 2)
    all_freq_df['donors'] = all_freq_df['0'] + all_freq_df['1'] + all_freq_df['2']
    all_freq_df['alle_freq'] = all_freq_df['sum'] / all_freq_df['donors']
    print(all_freq_df)
    # for col in ['NA', '0', '1', '2', 'sum', 'donors', 'alle_freq']:
    #     print(f"{col} - {set(all_freq_df[col])}")
    all_freq_df.to_csv(all_freq_path, sep="\t", index=False, encoding='utf-8', 
                compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})

    
    

def main():
    
    # File name vcf file
    name_vcf = "D:/Hanze_Groningen/STAGE/NEW PLAN/test_vcf.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/test_vcf.tsv.gz" # "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/" #all_vcf.tsv.gz" #"D:/Hanze_Groningen/STAGE/NEW PLAN/test_vcf.tsv.gz"
    all_freq_path = "D:/Hanze_Groningen/STAGE/NEW PLAN/all_freq.tsv.gz"
    df = pd.read_csv(name_vcf, sep='\t', compression='gzip')
    nan_value = False
    calculate_all_freq(df, all_freq_path, nan_value)

    
    


if __name__ == '__main__':
    main()




