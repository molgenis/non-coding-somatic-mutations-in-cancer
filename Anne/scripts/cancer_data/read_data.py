import pandas as pd


path_BRCA = 'C:/Users/Anne_/Downloads/simple_somatic_mutation.open.BOCA-UK.tsv'
df_BRCA = pd.read_csv(path_BRCA, sep='\t')

path_ORCA = 'C:/Users/Anne_/Downloads/simple_somatic_mutation.open.ORCA-IN.tsv'
df_ORCA = pd.read_csv(path_ORCA, sep='\t')

file_object = open('C:/Users/Anne_/Downloads/sample.tsv', 'a')
file_object.write(f'hoi\t{len(set(df_BRCA["icgc_donor_id"]))}')
file_object.close()



