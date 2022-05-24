#!/usr/bin/env python3
"""
Preprocessing to get tested and verified non-coding somatic snps

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
__title__ = "Preprocessing to get tested and verified non-coding somatic snps" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-24"
__updated__ = "2022-05-24"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 0.1
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse

import pandas as pd
from pybedtools import BedTool as bt


def main(args):
  df = pd.read_csv(args.inputPath, sep="\t")
  mutation_ids = pd.read_csv(args.mutationIDsPath, sep="\t", index_col=False, header=None)[0]
  print(mutation_ids)
  df = df[df['icgc_mutation_id'].isin(mutation_ids.tolist())]
  df = df[df['verification_status'] == "tested and verified"]
  non_coding_regions = ('intergenic_region', 'upstream_gene_variant', '5_prime_UTR_variant',
                        '3_prime_UTR_variant', 'downstream_gene_variant')
  df = df[df['consequence_type'].isin(non_coding_regions)]
  print(df)
  print(f"unique mutations = {df['icgc_mutation_id'].nunique()}")
  print(f"unique donors = {df['icgc_donor_id'].nunique()}")
  
  output_df = df[['chromosome','chromosome_start', 'chromosome_end', 'icgc_donor_id']]
  output_bt = bt.from_dataframe(output_df.rename(columns={'chromosome': '#chromosome'}), header=True)
  output_bt.saveas(args.outputPath)


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-m", "--mutationIDsPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Path to export the gene expression to")
  args = parser.parse_args()
  
  main(args)
