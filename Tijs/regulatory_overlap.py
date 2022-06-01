#!/usr/bin/env python3
"""
This script can be used to compare the overlap between SNPs and genome segmentation files in BED-format. 

MIT License

Copyright (c) 2022 Tijs van Lieshout

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Uses:
./regulatory_overlap.py -s 2022-04-22_eqtl_v1013_lead_snp_gene_with_info.bed -r 2022-04-13_GRCh37_UCNE.bed -o 2022-04-25_overlap_eqtl_snps_UCNE.bed --sortSNPs

./regulatory_overlap.py -s 2022-04-22_eqtl_v1013_lead_snp_gene_with_info.bed -r 2022-04-13_GRCh37_DNase.bed -o 2022-04-25_overlap_eqtl_snps_DNase.bed --sortSNPs

./regulatory_overlap.py -s 2022-04-22_eqtl_v1013_lead_snp_gene_with_info.bed -r 2022-04-13_GRCh37_TFBS.bed -o 2022-04-25_overlap_eqtl_snps_TFBS.bed --sortSNPs

./regulatory_overlap.py -s 2022-04-22_eqtl_v1013_lead_snp_gene_with_info.bed -r 2022-04-13_GRCh37_GREEN-DB.bed -o 2022-04-25_overlap_eqtl_snps_GREEN-DB.bed --sortSNPs


./regulatory_overlap.py -s ../ICGC_blood_data/tested_and_verified/2022-05-11_tested_and_verified_blood_snps.bed -r 2022-04-13_GRCh37_UCNE.bed.gz -o 2022-05-11_overlap_somatic_snps_UCNE.bed --sortSNPs

./regulatory_overlap.py -s ../ICGC_blood_data/tested_and_verified/2022-05-11_tested_and_verified_blood_snps.bed -r 2022-04-13_GRCh37_DNase.merged.bed.gz -o 2022-05-11_overlap_somatic_snps_DNase.bed --sortSNPs

./regulatory_overlap.py -s ../ICGC_blood_data/tested_and_verified/2022-05-11_tested_and_verified_blood_snps.bed -r 2022-04-13_GRCh37_TFBS.merged.bed.gz -o 2022-05-11_overlap_somatic_snps_TFBS.bed --sortSNP


./regulatory_overlap.py -s ../ICGC_blood_data/tested_and_verified/2022-05-11_tested_and_verified_blood_snps.bed -r 2022-04-13_GRCh37_DNase.merged.bed.gz -o 2022-05-11_overlap_somatic_snps_DNase.bed -c 2022-05-02_eqtl_v1013_lead_snp_gene_with_info.bed --sortSNPs


./regulatory_overlap.py -s ../../ICGC_blood_data/tested_and_verified/2022-05-12_non-coding_tested_and_verified_blood_snps.bed -r ../../chromatin_states/2022-05-12_healthy_blood_H3K4me1.bed -o 2022-05-13_overlap_somatic_snps_H3K4me1.bed -c ../../GREEN_DB/2022-05-02_eqtl_v1013_lead_snp_gene_with_info.bed --sortSNPs


./regulatory_overlap.py -s ../../GREEN_DB/2022-05-19_eqtl_v1013_lead_snp_gene_symbol_added.bed -r ../../GREEN_DB/2022-04-13_GRCh37_TFBS.merged.bed.gz -o 2022-05-30_overlap_stratified_eqtl_snps_TFBS.bed -e ../../GREEN_DB/2022-05-19_cis-eGenes_sorted_on_expr.tsv --sortSNPs
"""

# Metadata
__title__ = "Calculate the overlap between SNPs and regulatory regions" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-04-25"
__updated__ = "2022-06-01"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 1.4
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""
# Imports
import argparse

from pybedtools import BedTool as bt
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib_venn import venn2_unweighted, venn2_circles
from textwrap import wrap
from scipy.stats import chi2_contingency
import statsmodels.api as sm

from utilities import preprocess_bed_file


def main(args):
  snps_bed = preprocess_bed_file(args.SNPs) 
  reg_bed = preprocess_bed_file(args.RegRegion, True)
  if args.ComparisonSNPs:
    comparison_snps_bed = preprocess_bed_file(args.ComparisonSNPs) 

  if args.sortSNPs:
    snps_bed = snps_bed.sort()
    if args.ComparisonSNPs:
      comparison_snps_bed = comparison_snps_bed.sort()

  if args.sortRegRegion:
    reg_bed = reg_bed.sort()

  if args.ExpressionFile:
    snps_df = add_mean_gene_expression(snps_bed.to_dataframe(), args.ExpressionFile)
    
    contingency_df = pd.DataFrame()
    for q in np.sort(pd.unique(snps_df['quantile'])):
      subset_df = snps_df[snps_df['quantile'] == q]
      print(q)
      subset_bed = bt.from_dataframe(subset_df[['chrom', 'start', 'end']]).sort()
      overlap, jaccard = compute_jaccard(subset_bed, reg_bed)
      print('\n')

      n_overlap = jaccard['n_intersections']
      n_no_overlap = len(subset_bed) - jaccard['n_intersections']
      contingency_df[q] = [n_overlap, n_no_overlap]
    
    contingency_table = sm.stats.Table(contingency_df.values)
    CA_pvalue = contingency_table.test_ordinal_association().pvalue
    print(f"CA trend test p-value: {CA_pvalue:.2e}")
    return

  overlap, jaccard = compute_jaccard(snps_bed, reg_bed)
  compute_overlap(snps_bed, reg_bed, args.OutFile, args.RegRegion)

  if args.ComparisonSNPs:
    comparison_overlap, comparison_jaccard = compute_jaccard(comparison_snps_bed, reg_bed)
    p_value = calculate_chi2(jaccard['n_intersections'], comparison_jaccard['n_intersections'],
                             len(snps_bed), len(comparison_snps_bed))
    print(f"p-value = {p_value}")
    visualize_overlap_comparison(overlap, comparison_overlap, args.SNPs, 
                                 args.ComparisonSNPs, args.RegRegion, args.OutFile, p_value)


def add_mean_gene_expression(df, expr_path):
  expr_df = pd.read_csv(expr_path, sep='\t') 
  if 'score' in df.columns:
    column_name  = 'score'
  else:
    column_name = 'name'
  df = df[df[column_name].isin(expr_df['gene_symbol'])]
  merged_df = df.merge(expr_df, left_on=column_name, right_on='gene_symbol')

  # highest expressions are found in Quantile 1
  merged_df['quantile'] = pd.qcut(merged_df['mean_exp'], 4, labels=['Q4', 'Q3', 'Q2', 'Q1'])
  
  return merged_df


def compute_jaccard(snps_bed, reg_bed):
  jaccard = snps_bed.jaccard(reg_bed)
  overlap =(jaccard['n_intersections'] / len(snps_bed)) * 100

  print("Resulting overlap between SNPs and regulatory region:")
  for k, v in jaccard.items():
    print(f"{k:>15}{v:>15}")
  print()
  print(f"Which is {overlap:.2f}% of all SNPs")
  return overlap, jaccard


def compute_overlap(snps_bed, reg_bed, output_path, reg_path):
  overlap = snps_bed.intersect(reg_bed, wa=True, wb=True, filenames=True)
  overlap.saveas(output_path)

  venn2_unweighted(subsets = (len(snps_bed), len(reg_bed), len(overlap)),
                   set_labels = ("SNPs", reg_path.split("_")[-1]),
                   set_colors = ('orange', 'blue'),
                   alpha=0.7)

  plt.tight_layout()
  plt.savefig(f'{output_path}.png', dpi=300)
  plt.show()


def calculate_chi2(n_intersect, comparison_n_intersect, n_snps, comparison_n_snps):
  chi2, p_value, dof, expected = chi2_contingency(np.array([[n_intersect, comparison_n_intersect],
                                                            [n_snps, comparison_n_snps]]),
                                                  lambda_="log-likelihood")
  return p_value


def visualize_overlap_comparison(overlap, comparison_overlap, snps_path, comparison_snps_path, reg_path, output_path, p_value):
  fig, ax = plt.subplots(figsize=(4, 4))
  labels = [snps_path.split("/")[-1].split(".")[0].replace("_", " "), 
            comparison_snps_path.split("/")[-1].split(".")[0].replace("_", " ")]
  labels = ['\n'.join(wrap(l, 20)) for l in labels]
  rects = ax.bar(labels, [overlap, comparison_overlap], color='black')

  for rect in rects:
            height = rect.get_height()
            ax.annotate(xy=(rect.get_x() + rect.get_width()/2., height),
                        text=f"{round(height, 1)}%", ha='center', va='bottom')

  for spine in ax.spines:
    ax.spines[spine].set_visible(False)

  title = '\n'.join(wrap(f"proportion of SNPs in {reg_path.split('_')[-1].split('.')[0]}", 40))
  ax.set_title(title, pad=30)
  ax.annotate(xy=(0.5, 1.05), xycoords='axes fraction', text=f"p-value of chi2 < {p_value:.2e}",
              ha='center', va='bottom')
  ax.set(yticklabels=[])
  ax.tick_params(left=False)

  plt.tight_layout()
  output_title = title.split('\n')[-1]
  plt.savefig(f"{output_path}_comparison.png".replace(" ","_"), dpi=300)
  plt.show()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--SNPs", type=str, required=True, help="Path to BED of SNPs of interest")
  parser.add_argument("-r", "--RegRegion", type=str, required=True, help="Path to BED of regulatory regions of interest")
  parser.add_argument("-o", "--OutFile", type=str, required=True, help="Path to export the overlap BED of regulatory regions of interest with SNPs")
  parser.add_argument("--sortSNPs", action="store_true", required=False, help="Sort the BED of SNPs of interest LEXOGRAPHICALLY") 
  parser.add_argument("--sortRegRegion", action="store_true", required=False, help="Sort the BED of regulatory regions of interest LEXOGRAPHICALLY")
  parser.add_argument("-c", "--ComparisonSNPs", type=str, required=False, help="Path to BED of SNPs of interest to compare against")
  parser.add_argument("-e", "--ExpressionFile", type=str, required=False, help="Path to CSV of mean gene expression per eGene")
  args = parser.parse_args()
  
  main(args)
