#!/usr/bin/env python3
"""
Calculate the non-coding somatic mutational score per gene in multiple ways

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
./calculate_mutational_score_per_gene.py -i ../../ICGC_blood_data/2022-07-21_1200_upstream_of_TSS_blood_cancer_noncosomu_with_gene_symbols.bed -o ../../ICGC_blood_data/2022-07-22_stouffers_lymphoma_promoters.tsv
"""

# Metadata
__title__ = "Calculate the non-coding somatic mutational score per gene in multiple ways" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-07-22"
__updated__ = "2022-07-22"
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
import numpy as np
import scipy.stats as st


def main(args):
  df = pd.read_csv(args.inputPath, sep="\t")
  print(df)
  print("\n\n")


  # most naive way
  calculate_naive_count(df, 'icgc_mutation_id')
  print('positive strand SuRE scores')
  calculate_naive_score(df, 'SuRE_diff_positive')

  print('negative strand SuRE scores')
  calculate_naive_score(df, 'SuRE_diff_negative')

  # normalize mutational burden per donor
  normalized_count_matrix = create_normalized_matrix(df.groupby(['icgc_donor_id', 'name'])['icgc_mutation_id'].count())

  print('most naive way to sum Z-scores of positive strand SuRE p-values')
  normalized_pos_score_matrix = create_normalized_matrix(df.groupby(['icgc_donor_id', 'name'])['SuRE_diff_positive'].mean())
  print('most naive way to sum Z-scores of negative strand SuRE p-values')
  normalized_neg_score_matrix = create_normalized_matrix(df.groupby(['icgc_donor_id', 'name'])['SuRE_diff_negative'].mean())

  # calculate the stouffers methodology for combining z-scores
  calculate_stouffers(normalized_count_matrix)
  print('Stouffers Z-scores of positive strand SuRE p-values')
  calculate_stouffers(normalized_pos_score_matrix)
  print('Stouffers Z-scores of negative strand SuRE p-values')
  calculate_stouffers(normalized_neg_score_matrix)
#  stouffers_count.sort_values(ascending=False).to_csv(args.outputPath, sep="\t")


def calculate_naive_count(df, column_name):
  naive_count = df.groupby(['name'])[column_name].count()
  print("Most naive way to count mutational burden")
  print(naive_count.sort_values(ascending=False))
  print("\n\n")


def calculate_naive_score(df, column_name):
  naive_score = df.groupby(['name'])[column_name].sum()
  print("Most naive way to sum pred. scores")
  print(naive_score.sort_values(ascending=False))
  print("\n\n")


def create_normalized_matrix(per_donor):
  normalized_count_matrix = pd.DataFrame() 
  for donor in set(per_donor.index.get_level_values('icgc_donor_id')): 
    donor_data = per_donor[donor]
    donor_mean = donor_data.mean() 
    donor_std = donor_data.std()
    donor_norm = (donor_data - donor_mean) / donor_std
    donor_norm = donor_norm.rename(donor)
    normalized_count_matrix = pd.concat([normalized_count_matrix, donor_norm], axis=1)

  calculate_normalized_count(normalized_count_matrix)

  return normalized_count_matrix


def calculate_normalized_count(normalized_count_matrix):
  normalized_count = normalized_count_matrix.sum(axis=1)
  print("Sum of mutational counts normalized per patient (sum of Z-scores)")
  print(normalized_count.sort_values(ascending=False))
  print("\n\n")


def calculate_stouffers(normalized_matrix):
  # https://en.wikipedia.org/wiki/Fisher%27s_method#Relation_to_Stouffer.27s_Z-score_method
  k = len(normalized_matrix.columns) - normalized_matrix.isnull().sum(axis=1)
  stouffers_count = normalized_matrix.sum(axis=1) / np.sqrt(k) 
  print("Stouffer's Z-score method of combining mutational counts normalized per patient (Stouffer's method of Z-scores)")
  print(stouffers_count.sort_values(ascending=False))
  print("\n\n")


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
