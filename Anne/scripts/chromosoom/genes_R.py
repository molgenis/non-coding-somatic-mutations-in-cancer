import pandas as pd
import sys
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config

config = get_config()
gene_path = config['all_genes'] #'D:/Hanze_Groningen/STAGE/db/all_genes_new - kopie.tsv'
gene_df = pd.read_csv(gene_path, sep='\t')
print(gene_df.columns)
columns_select = ['#hg19.knownCanonical.chrom', 'hg19.knownCanonical.chromStart', 'hg19.knownCanonical.chromEnd', 'hg19.kgXref.geneSymbol'] #hg19.kgXref.geneSymbol   hg19.knownGene.name
select_gene_df = gene_df[columns_select].rename({'#hg19.knownCanonical.chrom': 'chr', 'hg19.knownCanonical.chromStart': 'start', 'hg19.knownCanonical.chromEnd': 'end', 'hg19.kgXref.geneSymbol': 'name'}, axis=1)
select_gene_df['gieStain'] = 'gpos50'

select_gene_df = select_gene_df[(~ select_gene_df['chr'].str.contains("_")) & (select_gene_df['chr'] != 'chrM')]
print(select_gene_df.head())
select_gene_df = select_gene_df[(~ select_gene_df['name'].str.contains("_"))]
print(set(select_gene_df['chr']))
# select_gene_df['name'] = select_gene_df['name'].str.replace('-','_')
select_gene_df.to_csv('D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt', sep='\t', encoding='utf-8', index = False)
