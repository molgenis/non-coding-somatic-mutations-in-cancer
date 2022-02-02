import os
import pandas as pd


file_name = 'donor' #'specimen' # 'donor'
your_path = "E:/STAGE/WGS/"
df = pd.DataFrame(columns=list(pd.read_csv(f'E:/STAGE/WGS/ALL-US/{file_name}.tsv', sep='\t').columns))
for path, dirs, files in os.walk(your_path):
    if f'{file_name}.tsv' in files:
        print(f'{path}/{file_name}.tsv')
        df_add = pd.read_csv(f'{path}/{file_name}.tsv', sep='\t') 
        df = pd.concat([df, df_add])
    else:
        print(f'---------------{path}')

df.to_csv(f'E:/STAGE/WGS/all_{file_name}.tsv', sep='\t', encoding='utf-8')

# print(list(pd.read_csv(f'E:/STAGE/WGS/ALL-US/{file_name}.tsv', sep='\t').columns))

# projects = list(set(df['icgc_donor_id']))
# print(projects)
# print(projects.sort())
# print(len(set(df['icgc_donor_id'])))