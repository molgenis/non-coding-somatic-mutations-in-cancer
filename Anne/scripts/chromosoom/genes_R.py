import pandas as pd


gene_path = 'D:/Hanze_Groningen/STAGE/db/all_genes_new - kopie.tsv' #'/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/cancer_data/genes_eQTL_etc/all_genes_new.tsv'
gene_df = pd.read_csv(gene_path, sep='\t')
select_gene_df = gene_df.iloc[:, 0:4].rename({'#hg19.knownCanonical.chrom': 'chr', 'hg19.knownCanonical.chromStart': 'start', 'hg19.knownCanonical.chromEnd': 'end', 'hg19.kgXref.geneSymbol': 'name'}, axis=1)
select_gene_df['gieStain'] = 'gpos50'
print(select_gene_df.head())
select_gene_df.to_csv('D:/Hanze_Groningen/STAGE/db/mycytobands_R.txt', sep='\t', encoding='utf-8', index = False)
