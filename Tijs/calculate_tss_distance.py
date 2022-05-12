#!/usr/bin/env python3
"""
Calculate the distance between TSS and SNPs

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
./calculate_tss_distance.py -s 2022-04-22_LEXOGRAPHICALLY_SORTED_eqtl_v1013_lead_snp_gene_with_info.bed.gz -t /groups/umcg-fg/tmp01/projects/non-coding-somatic/TSS_eQTL_distance/2022-05-09_TSS_hg19_canonical_ref.bed -o ../TSS_eQTL_distance/2022-05-09_TSS_eQTL_distance.bed

./calculate_tss_distance.py -s ../../ICGC_blood_data/tested_and_verified/2022-05-12_non-coding_tested_and_verified_blood_snps.bed -t ../../TSS_eQTL_distance/2022-05-09_TSS_hg19_canonical_ref.bed -o ../../TSS_eQTL_distance/2022-05-12_TSS_somatic_snps_distance.bed
"""

# Metadata
__title__ = "Calculate the distance between TSS and SNPs" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-04"
__updated__ = "2022-05-12"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 1.2
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse

from pybedtools import BedTool as bt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from utilities import preprocess_bed_file


def main(args):
  tss_df = sort_tss_bed(args.TSS)
  snps_bed = preprocess_bed_file(args.SNPs, True)
  tss_bed = bt.from_dataframe(tss_df, header=True)
  closest_df = compute_tss_distance(snps_bed, tss_bed, args.OutFile)
  plot_tss_distance(closest_df, args.OutFile)


def sort_tss_bed(tss_path):
  tss_df = pd.read_csv(tss_path, sep="\t", header=0, names=['#chrom', 'chromStart', 'chromEnd', 'name'],
                       dtype=str)
  tss_df['chromStart'] = pd.to_numeric(tss_df['chromStart'], errors='coerce')
  tss_df['chromEnd'] = pd.to_numeric(tss_df['chromEnd'], errors='coerce')
  tss_df = tss_df.dropna()
  tss_df['chromStart'] = tss_df['chromStart'].astype(int)
  tss_df['chromEnd'] = tss_df['chromEnd'].astype(int)
  tss_df = tss_df.sort_values(by=['#chrom', 'chromStart'])
  return tss_df


def compute_tss_distance(snps_bed, tss_bed, output_path):
  # report distance, keep only the first match between snp and tss
  closest_bed = tss_bed.closest(snps_bed, d=True, t="first")
  closest_bed.saveas(output_path)
  closest_df = closest_bed.to_dataframe()  
  #closest_df['manual_dist'] = abs(closest_df['thickStart'] - closest_df['start'])
  #print(closest_df.sample(30))

  print(100/len(closest_df) * len(closest_df[closest_df['itemRgb'] < 10000]))
  return closest_df
  

def plot_tss_distance(closest_df, output_path):
  plot_tss_distance_binned(closest_df, output_path)

  plot_tss_distance_whole(closest_df, output_path)


def plot_tss_distance_binned(closest_df, output_path):
  fig, ax = plt.subplots(figsize=(12, 6))

  bins = np.arange(0, 210000, 10000)
  ax.hist(np.clip(closest_df['itemRgb'], 0, bins[-1]), bins=bins, zorder=3)  
  xlabels = bins[1:].astype(str)
  xlabels = ["<" + str(x).replace("000", "") for x in xlabels] 
  xlabels[-1] = xlabels[-1].replace("<", ">")
  plt.xlim([0, 210000])
  plt.xticks(10000 * np.arange(len(xlabels)) + 5000)
  ax.set_xticklabels(xlabels) 
  ax.set_xlabel("Distance of SNP towards nearest gene (kb)", fontsize=18)
  ax.set_ylabel("No. genes", fontsize=18)

  ax.grid(color='lightgrey', axis='y', which='major', zorder=-3)
  ax.tick_params('x', labelbottom=True)
  for spine in ax.spines:
    ax.spines[spine].set_visible(False)

  plt.tight_layout()
  plt.savefig(f'{output_path}_BINNED.png', dpi=300)
  plt.show()


def plot_tss_distance_whole(closest_df, output_path):
  fig, ax = plt.subplots(figsize=(12, 6))

  n = len(closest_df['itemRgb'])
  x = np.linspace(closest_df['itemRgb'].min(), closest_df['itemRgb'].max())
  y = gaussian_kde(closest_df['itemRgb'])(x) * n
  ax.plot(x, y,
          color='black', zorder=3,)
  median = closest_df['itemRgb'].median()

  ax.vlines(median, 0, gaussian_kde(closest_df['itemRgb'])(median) * n,
            color='black', zorder=3)
  ax.annotate(xy=(median, gaussian_kde(closest_df['itemRgb'])(median) * n),
                  text=f"{round(median/1000, 0):.0f} kb", ha='left', va='bottom')
  ax.fill_between(x, y, alpha=0.25,
                  color='black', zorder=3)

  ax.set_xlabel("Distance of SNP towards nearest gene (kb)", fontsize=18)
  ax.set_ylabel("scaled density", fontsize=18)

  for spine in ax.spines:
    ax.spines[spine].set_visible(False)
  ax.grid(color='lightgrey', axis='y', which='major')
  ax.tick_params('x', labelbottom=True)

  plt.tight_layout()
  plt.savefig(f'{output_path}_WHOLE_median.png', dpi=300)
  plt.show()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--SNPs", type=str, required=True, help="Path to BED of SNPs of interest")
  parser.add_argument("-t", "--TSS", type=str, required=True, help="Path to BED of Transcription Start Sites of interest")
  parser.add_argument("-o", "--OutFile", type=str, required=True, help="Path to export the overlap BED of regulatory regions of interest with SNPs")
  args = parser.parse_args()
  
  main(args)
