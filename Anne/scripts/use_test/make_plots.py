import pandas as pd
import gzip
import io
import os

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

df = read_vcf("D:/Hanze_Groningen/STAGE/VCF/0001.vcf")
df_only_samples = df.iloc[:250,df.columns.get_loc("FORMAT"):df.columns.get_loc("FORMAT")+2]
print(df_only_samples.columns)

set_format = set()
for index, row in df_only_samples.iterrows():
    list_format = row['FORMAT'].split(':')
    set_format.update(list_format)
print(set_format)

df_format_value = pd.DataFrame(columns = [list(set_format)])
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
            list_values.append('.')
            # dict_row[i] = '.'
    df_format_value.loc[len(df_format_value)] = list_values
    # df_format_value = df_format_value.append(dict_row, ignore_index=True)

        
print(df_format_value.head())

