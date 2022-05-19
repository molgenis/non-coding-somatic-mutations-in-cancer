#!/usr/bin/env python3
"""
<A single line describing this program goes here.>

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
<The terminal interactions with this script go here>
"""

# Metadata
__title__ = "Template for a CLI python script" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-19"
__updated__ = "2022-05-19"
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


def main(args):
  exp_df = preprocess_exp(args.inputPath)
  exp_df.to_csv(args.outputPath, sep="\t", index=False)


def preprocess_exp(input_path):
  exp_df = pd.read_csv(input_path, sep=",")
  cis_exp_df = exp_df[exp_df['cis_eGene'] == 'yes'].sort_values(by=['mean_exp'], ascending=False)
  gene_symbol_overlap = len(cis_exp_df['gene_symbol'].dropna()) / len(cis_exp_df)
  print(f"overlap of ENSG with gene symbols = {gene_symbol_overlap}")
  cis_exp_df = cis_exp_df[['gene_symbol', 'mean_exp']].dropna()
  print(cis_exp_df)
  return cis_exp_df


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Path to export the gene expression to")
  args = parser.parse_args()
  
  main(args)
