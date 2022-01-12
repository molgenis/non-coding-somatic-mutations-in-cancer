import pandas as pd

# snp132_ucsc_hg19_UCSC, snp132_ucsc_hg19_NCBI, snp132_ucsc_hg19_V38
filename= 'snp132_ucsc_hg19_V38' 
df = pd.read_csv(f'D:/Hanze_Groningen/STAGE/bed/{filename}.bed', sep='\t')
print(df.head())
df.to_csv(f"D:/Hanze_Groningen/STAGE/bed/{filename}_gene.bed", sep="\t", index=False)
df['exonStarts'] = df['exonStarts'].str[:-1].str.strip().str.split(',')
df['exonEnds'] = df['exonEnds'].str[:-1].str.strip().str.split(',')
print(df.head())
#https://stackoverflow.com/questions/12680754/split-explode-pandas-dataframe-string-entry-to-separate-rows
#for multiple columns (for Pandas 1.3.0+):
df2 = df.explode(['exonStarts', 'exonEnds'])
print(df2.head())

df2.to_csv(f"D:/Hanze_Groningen/STAGE/bed/{filename}_exon.bed", sep="\t", index=False)