import pandas as pd


def compare_eQTL_files(df_eQTL, df_strong_eQTL):
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
    df.to_csv("/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/genes_eQTL_etc/eqtl_v1013_lead_snp_gene_with_info.txt", sep='\t', encoding='utf-8', index=False)




def main():
    path_eQTL = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/genes_eQTL_etc/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt"
    df_eQTL = pd.read_csv(path_eQTL, sep='\t')
    print('read 1 done')
    path_strong_eQTL = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/eQTL/eqtl_v1013_lead_snp_gene.txt"
    df_strong_eQTL = pd.read_csv(path_strong_eQTL, sep='\t', header=None)
    print('read 2 done')
    compare_eQTL_files(df_eQTL, df_strong_eQTL)



if __name__ == '__main__':
    main()