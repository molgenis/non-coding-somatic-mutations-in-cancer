#!/usr/bin/env python3
from numpy.core.numeric import NaN
import pandas as pd
import io
import sys
import os


# https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
def read_vcf(path):
    """

    :param path:
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

path = sys.argv[1].strip()
if ".gz" in path:
    path_command_file = path
    path = os.path.splitext(path)[0]
    num_samples = 2
else:
    path_command_file = path + '.gz'
    num_samples = 1
    
# Read vcf file
df = read_vcf(path)#(sys.argv[1].strip())
df_columns = list(df.columns[:-int(num_samples)])
name_new_file = '_'.join(list(df.columns[-int(num_samples):]))


# Create a set of all abbreviations in the FORMAT column that are separated by :
list_format = [i.split(':') for i in list(set(list(df['FORMAT'])))]
# Use list comprehension to convert a list of lists to a flat list 
set_format = set([ item for elem in list_format for item in elem])
# Create a set of all abbreviations in the INFO column that are separated by ;
# Replace all values =NUMBER; with ;
df['INFO'] = df['INFO'].replace(to_replace ='=(.*?)\;', value = ';', regex = True)
# Replace values =NUMBER with ''
set_info_joined = set(list(df['INFO'].replace(to_replace ='=.*', value = '', regex = True)))
list_info = [i.split(';') for i in list(set_info_joined)]
set_info = set([ item for elem in list_info for item in elem])
#print(set_format)
#print(set_info)

#https://www.biostars.org/p/226965/
#'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AF\n'
#bcftools query -Hf "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FORMAT\t%FILTER$FIELDS[\t%GT]\n" V42.vcf.gz > RESULTSSSexcel
stringToGetFiles= ''
for col_name in df_columns:
    if col_name not in ['INFO', 'FORMAT']:
        stringToGetFiles+=f'%{col_name}\\t'
for col_name in sorted(list(set_info), key=str.lower):
    stringToGetFiles+=f'%INFO/{col_name}\\t'
stringToGetFiles = stringToGetFiles[:-2]
for col_name in sorted(list(set_format), key=str.lower):
    stringToGetFiles+=f'[\\t%{col_name}]' #[\t%GT]
stringToGetFiles+='\\n'

# bcftools query -Hf "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/CONTQ\t%INFO/DP\t%INFO/ECNT\t%INFO/GERMQ\t%INFO/MBQ\t%INFO/MFRL\t%INFO/MMQ\t%INFO/MPOS\t%INFO/POPAF\t%INFO/RPA\t%INFO/RU\t%INFO/SEQQ\t%INFO/STR\t%INFO/STRANDQ\t%INFO/STRQ\t%INFO/TLOD[\t%AD][\t%AF][\t%DP][\t%F1R2][\t%F2R1][\t%GT][\t%PGT][\t%PID][\t%PS][\t%SB]\n" /groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/datasets/EGAD00001000292/samples/S6/compare_5042_5044/bowtie/0001.vcf.gz > 00011RESULTSSSexcel.vcf
part_command = f'"{stringToGetFiles}" {path_command_file} -o {sys.argv[2]}{name_new_file}.vcf'
#print(part_command)

f = open(sys.argv[3].strip(), 'a')
f.write(f'{part_command}\n')
f.close()