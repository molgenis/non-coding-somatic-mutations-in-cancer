#!/usr/bin/env python3
import pandas as pd
import io
import seaborn as sns
import matplotlib.pyplot as plt
import sys

#sys.argv[1] = path file
#sys.argv[2] = 
#sys.argv[3] = 
path = "/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S1/compare_4094_4099/bowtie/0001.vcf"

# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

df = read_vcf(sys.argv[1].strip())
df_only_samples = df.iloc[:,df.columns.get_loc("FORMAT"):]
name = df_only_samples.columns[-1]
print(name)
print(df_only_samples.columns)

set_format = set()
for index, row in df_only_samples.iterrows():
    list_format = row['FORMAT'].split(':')
    set_format.update(list_format)
print(set_format)

df_format_value = pd.DataFrame(columns = list(set_format))
for index, row in df_only_samples.iterrows():
    list_values=list()
    # dict_row = dict()
    names_format = row['FORMAT'].split(':')
    values_format = row[1].split(':')    
    for i in list(set_format):
        if i in names_format:
            index_value = names_format.index(i)
            # dict_row[i]=  values_format[index_value]  
            list_values.append(values_format[index_value] )
        else:
            list_values.append(float("NaN"))
            # dict_row[i] = '.'
    df_format_value.loc[len(df_format_value)] = list_values
    # df_format_value = df_format_value.append(dict_row, ignore_index=True)
print(df_format_value.head())
df_format_value.to_csv(f'{sys.argv[2]}_{name}.csv', sep='\t', encoding='utf-8')

# for column in list(df_format_value.columns):
#     data = df_format_value[df_format_value[column].notna()]
#     if column in ['AF', 'DP']:
#         print(data[column][25:30])
#         data[column] = pd.to_numeric(data[column], downcast="float")
#         if not data.empty:    
#             plt.figure(figsize=(15, 15))
#             sns.histplot(data=data, x=column).set_title(f'{column}')
#             plt.yscale('log')
#             plt.savefig(f'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/FORMAT/{column}_seaborn_hist.png')
#             plt.clf()
#         else:
#             print(f'{column} EMPTY')
