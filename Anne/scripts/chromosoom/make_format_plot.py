#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

# Also takes the folder 1 higher, so that I can do the import after
# sys.path.append("..")
sys.path.append('/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from vcf_reader import read_vcf


def make_plot_format_vcf(path, basename):
    basename = basename.split('.')[0]
    # Read vcf file
    df = read_vcf(path)#(sys.argv[1].strip())    
    samples_df = list(df.columns)[9:]
    general_df = df.iloc[:, 0:9]
    for sample in samples_df:
        general_df['sample'] = df[sample]
        filter_sample = general_df[~general_df['sample'].str.contains(":.:.")]
        filter_sample["numCHROM"] = pd.to_numeric(filter_sample["CHROM"].str.replace("chr", ""))
        sort_df = filter_sample.sort_values(["numCHROM", "POS"])
        sort_df['SNPnum'] = np.arange(len(sort_df))
        plot_format = sort_df[['SNPnum', 'CHROM', 'POS', 'POS']]
        plot_format.to_csv(f'D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/vcf/{sample}_{basename}.bed', sep="\t", index=False, header=None)

def make_plot_format_other(path, basename):
    basename = 'download?fn=%2Fcurrent%2FProjects%2FBOCA-UK%2Fsimple_somatic_mutation.open.BOCA-UK.tsv'
    basename = basename.split('%2F')[3]
    df = pd.read_csv(path, sep='\t')
    samples_df = list(set(df['icgc_donor_id']))
    print(len(samples_df))
    for sample in samples_df:
        print(f'----{sample}')
        filter_sample = df.loc[df['icgc_donor_id'] == sample]
        filter_sample = filter_sample.loc[filter_sample['sequencing_strategy'] == 'WGS']
        # Remove duplicates (de kolommen in de list kunnen wel anders zijn, maar worden dan toch verwijderd (eerste voorbeeld wordt dan gehouden))
        filter_sample = filter_sample.drop_duplicates(subset=filter_sample.columns.difference(['consequence_type', 'aa_mutation', 'cds_mutation', 'gene_affected', 'transcript_affected']))
        print(len(filter_sample))
        if len(filter_sample) > 10:            
            sort_df = filter_sample.sort_values(["chromosome", "chromosome_start"])            
            plot_format = sort_df[['chromosome', 'chromosome_start', 'chromosome_end']]
            #plot_format = plot_format.drop_duplicates() 
            plot_format.insert(loc=0, column='SNPnum', value=np.arange(len(plot_format)))
            print(len(plot_format))
            #sort_df['SNPnum'] = np.arange(len(sort_df))
            plot_format['chromosome'] = 'chr' + plot_format['chromosome'].astype(str)
            plot_format.to_csv(f'D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/other/{sample}_{basename}.bed', sep="\t", index=False, header=None)
            
        
path = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/simple_somatic_mutation.open.BOCA-UK.tsv" 
# path = "D:/Hanze_Groningen/STAGE/DIFFERENT CANCERS/merge_manual_bwa_aln.vcf"
print(path)
# Get the basename of the file
basename = os.path.basename(path) #.split('.')[0]
type_file = 'xxx'
if type_file == 'vcf':
    make_plot_format_vcf(path, basename)
else:
    make_plot_format_other(path, basename)