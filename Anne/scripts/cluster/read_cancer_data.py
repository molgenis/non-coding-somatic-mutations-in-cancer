import pandas as pd
import sys
import os

path = sys.argv[1] #.replace('.gz', '') #'C:/Users/Anne_/Downloads/simple_somatic_mutation.open.BRCA-US.tsv'

df = pd.read_csv(path, sep='\t')


lijst_col = list(df.columns)
#lijst_col.remove('quality_score')
#lijst_col.remove('probability')
#lijst_col.remove('verification_platform')
#lijst_col.remove('biological_validation_platform')
#lijst_col.remove('initial_data_release_date')
#lijst_col.remove('experimental_protocol')

print(path.split("%2F")[3])
print(set(df["sequencing_strategy"]))
print(len(set(df["icgc_donor_id"])))
print(len(df))

file_object = open(f'{sys.argv[2]}cancer_number5.tsv', 'a')

#if os.stat(f'{sys.argv[2]}cancer_number3.tsv').st_size == 0:
#    file_object.write(f'cancer  icgc_donor_id') #{"  ".join(list(df.columns[1]))}\n')
    
file_object.write(f'{path.split("%2F")[3]}')




#for col in list(df.columns[1]):
    #print(col)
    #print(len(set(df[col])))
file_object.write(f'\t{len(set(df["icgc_donor_id"]))}')
file_object.write(f'\t{set(df["sequencing_strategy"])}') 
file_object.write(f'\t{len(df)}') 

df = df.loc[df["sequencing_strategy"] == 'WGS']

file_object.write(f'\t{len(set(df["icgc_donor_id"]))}')
file_object.write(f'\t{set(df["sequencing_strategy"])}') 
print(set(df["sequencing_strategy"]))
print(len(set(df["icgc_donor_id"])))
print(len(df))
file_object.write(f'\t{len(df)}')    
file_object.write(f'\n')
file_object.close()