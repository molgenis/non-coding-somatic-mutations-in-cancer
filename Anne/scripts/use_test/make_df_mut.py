#!/usr/bin/env python3
import pandas as pd
import io
import sys

# path_file = "D:/Hanze_Groningen/STAGE/VCF/S1_numT_1_numHC_1_SS6004094_SS6004099__somatic_filtered.vcf"
# path_file = "D:/Hanze_Groningen/STAGE/VCF/0001.vcf"

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
hc_tum = df.columns[-2:]

df_only_samples = df.iloc[:200,df.columns.get_loc("FORMAT"):]

set_format = set()
for index, row in df_only_samples.iterrows():
    list_format = row['FORMAT'].split(':')
    set_format.update(list_format)
print(set_format)

df_format_value_hc = pd.DataFrame(columns = list(set_format))
df_format_value_tum = pd.DataFrame(columns = list(set_format))
for index, row in df_only_samples.iterrows():
    list_values_hc=list()
    list_values_tum=list()
    names_format = row['FORMAT'].split(':')
    values_hc_format = row[hc_tum[0]].split(':')
    values_tum_format = row[hc_tum[1]].split(':')     
    for i in list(set_format):
        if i in names_format:
            index_value = names_format.index(i)
            list_values_hc.append(values_hc_format[index_value])
            list_values_tum.append(values_tum_format[index_value])
        else:
            list_values_hc.append(float("NaN"))
            list_values_tum.append(float("NaN"))
    df_format_value_hc.loc[len(df_format_value_hc)] = list_values_hc
    df_format_value_tum.loc[len(df_format_value_tum)] = list_values_tum
    
df_format_value_hc.to_csv(f'{sys.argv[2]}{hc_tum[0].split(".")[0]}_{hc_tum[1].split(".")[0].split("_")[1]}_hc.csv', sep='\t', encoding='utf-8')
df_format_value_tum.to_csv(f'{sys.argv[2]}{hc_tum[1].split(".")[0]}.csv', sep='\t', encoding='utf-8')


