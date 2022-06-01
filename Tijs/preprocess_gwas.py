#!/usr/bin/env python3
"""
Preprocessing to get gwas non-coding somatic snps

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
./preprocess_icgc.py -i ../../ICGC_blood_data/tested_and_verified/2022-05-24_simple_somatic_mutation.open.tsv -m ../../ICGC_blood_data/tested_and_verified/2022-05-24_mutation_ids_of_interest.tsv -o ../../ICGC_blood_data/tested_and_verified/2022-05-24_tested_and_verified_non_coding_icgc.bed
"""

# Metadata
__title__ = "Preprocessing to get gwas non-coding snps" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-06-01"
__updated__ = "2022-06-01"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 1.0
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse

import pandas as pd
from pybedtools import BedTool as bt
import numpy as np


def main(args):
  df = pd.read_csv(args.inputPath, sep="\t")
  if args.dropCoding: 
    coding_regions = ('non-coding_transcript_exon_variant', 'missense_variant', 'synonymous_variant',
                      'stop_gained', 'frameshift_variant', 'start_lost', 'protein_altering_variant')
    df = df[~df['functional_class'].isin(coding_regions)]
  df = df[['chromosome_name', 'chromosome_position']].dropna()
  df['chromosome_position'] = df['chromosome_position'].astype(int)
  df['chromEnd'] = df['chromosome_position'] + 1
  df[df.columns[0]] = np.where(df[df.columns[0]].str.contains("chr"), 
                               df[df.columns[0]], 
                               'chr' + df[df.columns[0]])
  output_bt = bt.from_dataframe(df.rename(columns={'chromosome_name': '#chrom', 'chromosome_position': 'chromStart',}),
                                          header=True)
  output_bt.sort().saveas(args.outputPath)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Path to export the gene expression to")
  parser.add_argument("-d", "--dropCoding", action="store_true", required=False, help="Drop coding regions") 
  args = parser.parse_args()
  
  main(args)
