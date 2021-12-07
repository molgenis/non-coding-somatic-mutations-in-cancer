import pandas as pd
import sys
import os

path = sys.argv[1] #.replace('.gz', '') #'C:/Users/Anne_/Downloads/simple_somatic_mutation.open.BRCA-US.tsv'

df = pd.read_csv(path, sep='\t')

file_object = open(f'{sys.argv[2]}cancer_number2.tsv', 'a')

if os.stat(path).st_size == 0:
    file_object.write(f'cancer  {"  ".join(list(df.columns))}\n')
    
file_object.write(f'{path.split("%2F")[3]}')
for col in df.columns:
    file_object.write(f'\t{len(set(df[col]))}')
    
file_object.write(f'\n')
file_object.close()