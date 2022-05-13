import pandas as pd 
import ast


def prep_file(path_file, path_save, type_bef_aft):
    df = pd.read_csv(path_file, sep='\t')

    breast_count_list = list()
    nonbreast_count_list = list()
    for index, row in df.iterrows():
        print(index)
        if row['cancer_count'] != '-':
            dict_cancer_count = ast.literal_eval(row['cancer_count'])
            total_count = sum(dict_cancer_count.values())
            if 'Breast' in dict_cancer_count:
                breast_count_list.append(dict_cancer_count['Breast'])
                nonbreast_count = total_count - dict_cancer_count['Breast']
                nonbreast_count_list.append(nonbreast_count)
            else:
                breast_count_list.append(0)
                nonbreast_count_list.append(total_count)
        else:
            breast_count_list.append(0)
            nonbreast_count_list.append(0)
            
    df_b_nb = df.iloc[:, 1:5]
    df_b_nb['counts_breast'] = breast_count_list
    df_b_nb['counts_nonbreast'] = nonbreast_count_list
    print(df_b_nb)
    df_b_nb.to_csv(f"{path_save}b_nb_gene_{type_bef_aft}_2000_250.tsv", sep='\t', encoding='utf-8', index=False)
    return df_b_nb

def main():
    type_bef_aft = 'after'
    path_save = 'D:/Hanze_Groningen/STAGE/UMAP/'
    path_file = f"D:/Hanze_Groningen/STAGE/UMAP/ALL_gene_{type_bef_aft}_2000_250.tsv"
    df_b_nb = prep_file(path_file, path_save, type_bef_aft)



if __name__ == '__main__':
    main()

