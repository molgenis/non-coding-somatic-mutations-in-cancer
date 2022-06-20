import pandas as pd




def main():
    path_file_gene = "D:/Hanze_Groningen/STAGE/analyse/new/correction/run_cluster/DNase_NonCoding_all_MTC_close_gene.tsv"
    path_file_snps = "D:/Hanze_Groningen/STAGE/analyse/new/correction/run_cluster/DNase_NonCoding_sig_snps.tsv"
    df_gene = pd.read_csv(path_file_gene, sep='\t')
    print(df_gene)
    df_snps = pd.read_csv(path_file_snps, sep='\t')
    print(df_snps)






if __name__ == '__main__':
    main()