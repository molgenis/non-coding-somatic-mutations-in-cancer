#!/usr/bin/env python3

#Imports
import pandas as pd
import numpy as np
import sys

sys.path.append(
    '/groups/umcg-wijmenga/tmp01/projects/lude_vici_2021/rawdata/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def search_gene(list_chr, path_analyse, type_file, non_coding, type_MTC, fc):
    """
    Finds the closest gene to a significant region
    Also checks whether the region is in the gene, overlaps with the gene or whether the gene is completely
    outside the region
    :param list_chr: List with chromosomes
    :param path_analyse:  The path where the new files will be saved
    :param type_file:  Name of the layer 
    :param non_coding: Noncoding data, coding data or mix
    :param type_MTC: Type of multiple testing correction
    :param fc: Type fold change (all, fc<1 or fc> 1)  
    :return:    
    """
    if len(list_chr) > 0:
        # Call get_config
        config = get_config('gearshift')
        path_gene_file = config['all_genes']
        df = pd.read_csv(path_gene_file, sep='\t')
        # Open file
        f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{type_MTC}_{fc}_close_gene.tsv", "a")
        # Write header
        f.write(
            f"snp_ID\tchr\tpos_begin\tpos_end\tfc\tbigger\thg19.kgXref.geneSymbol\tposition_region_to_gene\t"
            f"distance_gene\tcolumn_name_min\n")
        for pos in list_chr:
            if len(pos.split('_')) == 5:
                chr, pos_begin, pos_end, fc, bigger = pos.split('_')
                snp_ID = '-'
            elif len(pos.split('__')) == 6:
                snp_ID, chr, pos_begin, pos_end, fc, bigger = pos.split('__')
            # Check if chr starts with 'chr'
            if chr.startswith("chr"):
                df_select = df.loc[df['hg19.knownGene.chrom'] == chr].reset_index(drop=True)
            else:
                chr = f'chr{chr}'
                df_select = df.loc[df['hg19.knownGene.chrom'] == chr].reset_index(drop=True)
            # in gene
            df_filter = df_select[(df_select['hg19.knownGene.txStart'] <= round(float(pos_begin))) & (
                        df_select['hg19.knownGene.txEnd'] >= round(float(pos_end)))]
            if len(df_filter) < 1:
                # overlap with gene
                df_filter = df_select[((df_select['hg19.knownGene.txStart'] >= round(float(pos_begin))) & (
                            df_select['hg19.knownGene.txStart'] <= round(float(pos_end))) & (
                                                   df_select['hg19.knownGene.txEnd'] >= round(float(pos_end)))) |
                                      ((df_select['hg19.knownGene.txStart'] <= round(float(pos_begin))) & (
                                                  df_select['hg19.knownGene.txEnd'] <= round(float(pos_end))) & (
                                                   df_select['hg19.knownGene.txEnd'] >= round(float(pos_begin)))) |
                                      ((df_select['hg19.knownGene.txStart'] >= round(float(pos_begin))) & (
                                                  df_select['hg19.knownGene.txEnd'] <= round(float(pos_end))))]
                if len(df_filter) < 1:
                    # out gene
                    df_filter = df_select
                    df_filter['num_txStart_begin'] = (
                                df_filter['hg19.knownGene.txStart'] - round(float(pos_begin))).abs()
                    df_filter['num_txStart_end'] = (df_filter['hg19.knownGene.txStart'] - round(float(pos_end))).abs()
                    df_filter['num_txEnd_begin'] = (df_filter['hg19.knownGene.txEnd'] - round(float(pos_begin))).abs()
                    df_filter['num_txEnd_end'] = (df_filter['hg19.knownGene.txEnd'] - round(float(pos_end))).abs()

                    list_min = [df_filter['num_txStart_begin'].loc[df_filter['num_txStart_begin'].idxmin()],
                                df_filter['num_txStart_end'].loc[df_filter['num_txStart_end'].idxmin()],
                                df_filter['num_txEnd_begin'].loc[df_filter['num_txEnd_begin'].idxmin()],
                                df_filter['num_txEnd_end'].loc[df_filter['num_txEnd_end'].idxmin()]]
                    colnames_min = ['num_txStart_begin', 'num_txStart_end', 'num_txEnd_begin', 'num_txEnd_end']
                    value_min = min(list_min)
                    ind = np.argmin(list_min)
                    index_min = df_filter[colnames_min[ind]].idxmin()
                    f.write(
                        f"{snp_ID}\t{chr}\t{pos_begin}\t{pos_end}\t{fc}\t{bigger}"
                        f"\t{df_filter.iloc[index_min]['hg19.kgXref.geneSymbol']}\tout\t{value_min}\t{colnames_min[ind]}\n")
                else:
                    f.write(
                        f"{snp_ID}\t{chr}\t{pos_begin}\t{pos_end}\t{fc}\t{bigger}\t"
                        f"{','.join(map(str, list(set(list(df_filter['hg19.kgXref.geneSymbol'])))))}\toverlap\t0\t-\n")
            else:
                f.write(
                    f"{snp_ID}\t{chr}\t{pos_begin}\t{pos_end}\t{fc}\t{bigger}\t"
                    f"{','.join(map(str, list(set(list(df_filter['hg19.kgXref.geneSymbol'])))))}\tin\t0\t-\n")
        f.close()
    else:
        f = open(f"{path_analyse}correction/{type_file}_{non_coding}_{type_MTC}_{fc}_close_gene.tsv", "a")
        f.write(
            f"snp_ID\tchr\tpos_begin\tpos_end\tfc\tbigger\thg19.kgXref.geneSymbol\t"
            f"position_region_to_gene\tdistance_gene\tcolumn_name_min\n")
        f.write(f'-\t-\t-\t-\t-\t-\t-\t-\t-\n')
        f.close()
