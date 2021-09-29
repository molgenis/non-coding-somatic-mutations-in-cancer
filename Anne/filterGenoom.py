"""

"""

import pandas as pd
import pysam

# read file
df = pd.read_table('Table_browser_19', sep='\t')
# Select columns for bed file
bed_file = df[['chrom', 'txStart', 'txEnd']]
# unique combinations of values in selected columns
uni_bed_file = bed_file.groupby(['chrom', 'txStart', 'txEnd']).size().reset_index().iloc[:, :-1]
# .rename(columns={0:'count'}).drop('count', axis=1)
# write bed file
# uni_bed_file.to_csv('test.bed', sep='\t', encoding='utf-8')
# for index, row in df.iterrows():
