import pandas as pd 

def sort_on_column_layers(df, column_name, top_num):
    sort_df = df.sort_values(by=[column_name])
    print('_______________________________________________')
    select_top = sort_df.head(top_num)
    print(select_top)
    return select_top
    
def print_statements(self_X2_df, X2_df, F_df, column_name):
    print('---ALL')
    all_values = set(self_X2_df[column_name]) | set(X2_df[column_name]) |  set(F_df[column_name])
    print(f'All genes ({len(all_values)}): {all_values}')
    elements_in_all = list(set.intersection(*map(set, [self_X2_df[column_name], X2_df[column_name], F_df[column_name]])))
    print(f'Gene in all list ({len(elements_in_all)}): {elements_in_all}\n')
    # selfX2 and X2
    print('---selfX2 and X2')
    in_selfX2_X2 = set(self_X2_df[column_name]) & set(X2_df[column_name])
    print(f'in selfX2 and X2 ({len(in_selfX2_X2)}): {in_selfX2_X2}')
    dif_selfX2_X2 = set(self_X2_df[column_name]) - set(X2_df[column_name])
    print(f'in selfX2 NOT in X2 ({len(dif_selfX2_X2)}): {dif_selfX2_X2}')
    dif_X2_selfX2 = set(X2_df[column_name]) - set(self_X2_df[column_name])
    print(f'in X2 NOT in selfx2 ({len(dif_X2_selfX2)}): {dif_X2_selfX2}\n')
    # selfX2 and F
    print('---selfX2 and F')
    in_selfX2_F = set(self_X2_df[column_name]) & set(F_df[column_name])
    print(f'in selfX2 and F ({len(in_selfX2_F)}): {in_selfX2_F}')
    dif_selfX2_F = set(self_X2_df[column_name]) - set(F_df[column_name])
    print(f'in selfX2 NOT in F ({len(dif_selfX2_F)}): {dif_selfX2_F}')
    dif_F_selfX2 = set(F_df[column_name]) - set(self_X2_df[column_name])
    print(f'in F NOT in selfX2 ({len(dif_F_selfX2)}): {dif_F_selfX2}\n')
    # F and X2
    print('---F and X2')
    in_F_X2 = set(F_df[column_name]) & set(X2_df[column_name])
    print(f'in F and X2 ({len(in_F_X2)}): {in_F_X2}')
    dif_F_X2 = set(F_df[column_name]) - set(X2_df[column_name])
    print(f'in F NOT in X2 ({len(dif_F_X2)}): {dif_F_X2}')
    dif_X2_F = set(X2_df[column_name]) - set(F_df[column_name])
    print(f'in X2 NOT in F ({len(dif_X2_F)}): {dif_X2_F}')

def layers(path_file, top_num, column_name):
    df = pd.read_csv(path_file, sep='\t')
    if column_name == 'snp':
        df['snp'] = df['chr'].map(str) + '_' + df['pos_start'].map(str) + '_' + df['pos_end'].map(str)
    select_top_X2_self = sort_on_column_layers(df, 'p_value_X2_self', top_num)
    select_top_X2 = sort_on_column_layers(df, 'p_value_X2', top_num)
    select_top_F = sort_on_column_layers(df, 'p_value_F', top_num)
    # print_statements(select_top_X2_self, select_top_X2, select_top_F, column_name)



def main():
    # # per_snp
    # path_file = "D:/Hanze_Groningen/STAGE/analyse/stat/per_snp_ALL_both_0_TESTS.tsv"
    # layers(path_file, 20, 'snp')
    # path_file = "D:/Hanze_Groningen/STAGE/analyse/stat/per_snp_Coding_both_0_TESTS.tsv"
    # layers(path_file, 20, 'snp')
    # path_file = "D:/Hanze_Groningen/STAGE/analyse/stat/per_snp_NonCoding_both_0_TESTS.tsv"
    # layers(path_file, 20, 'snp')


    # bef
    path_file = "D:/Hanze_Groningen/STAGE/analyse/stat/beforeGene_NonCoding_Coding_both_0_TESTS.tsv"
    layers(path_file, 20, 'gene')
    # aft
    # path_file = "D:/Hanze_Groningen/STAGE/analyse/stat/afterGene_NonCoding_Coding_both_0_TESTS.tsv"
    # layers(path_file, 20, 'gene')
    # DNase


    # TFBS


    # UCNE


    # GT
    

    
    


if __name__ == '__main__':
    main()