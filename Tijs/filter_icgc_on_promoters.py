#!/usr/bin/env python3
"""
Creating a set of all given ICGC mutations 1200 bases upstream of the canonical transcription start sites

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
./filter_icgc_on_promoters.py -i ../../ICGC_blood_data/2022-07-01_simple_somatic_mutation.open.tsv.gz -t ../../ICGC_blood_data/2022-05-09_TSS_hg19_canonical_ref.bed -o ../../ICGC_blood_data/2022-07-13_1200_upstream_of_TSS_blood_cancer_noncosomu.bed
"""

# Metadata
__title__ = "Creating a set of all given ICGC mutations 1200 bases upstream of the canonical transcription start sites"
__author__ = "Tijs van Lieshout"
__created__ = "2022-07-01"
__updated__ = "2022-07-13"
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

from tqdm import tqdm


def main(args):
  tss_df = sort_tss_bed(args.TSS)

  # file is too large to load in at once, requires loading in chunks
  chunk_size = 10 ** 6
  upstream_df = pd.DataFrame()
  for chunk in tqdm(pd.read_csv(args.inputPath, sep="\t", compression='gzip', chunksize=chunk_size,
                    dtype=str, usecols=['chromosome','chromosome_start', 'chromosome_end', 'icgc_donor_id', 
                    'icgc_mutation_id', 'mutated_from_allele', 'mutated_to_allele'])):
    # dtype=str causes everything to be string, format the required columns back to int:
    chunk['chromosome_start'] = chunk['chromosome_start'].astype(int)

    # do the main processing step
    processed_chunk = get_upstream_mutations(chunk, tss_df, 1200)
    upstream_df = pd.concat([upstream_df, processed_chunk])

  upstream_df = upstream_df.drop_duplicates()
  print(f"unique mutations = {upstream_df['icgc_mutation_id'].nunique()}")
  print(f"unique donors = {upstream_df['icgc_donor_id'].nunique()}")
  
  # prepare output, change column order
  output_df = upstream_df[['chromosome','chromosome_start', 'chromosome_end', 'icgc_donor_id', 
                           'icgc_mutation_id', 'mutated_from_allele', 'mutated_to_allele']]
  output_df = output_df.sort_values(by=['chromosome', 'chromosome_start'])

  output_bt = bt.from_dataframe(output_df.rename(columns={'chromosome': '#chromosome'}), header=True)
  output_bt.saveas(args.outputPath)


def sort_tss_bed(tss_path):
  """
  sort and preprocess bed file containing transcription start sites
  :param tss_path: path to bed file including the .bed part
  :return: dataframe equivalent of bed file, sorted on chromosome position
  """
  tss_df = pd.read_csv(tss_path, sep="\t", header=0, names=['#chrom', 'chromStart', 'chromEnd', 'name'],
                       dtype=str)
  tss_df['chromStart'] = pd.to_numeric(tss_df['chromStart'], errors='coerce')
  tss_df['chromEnd'] = pd.to_numeric(tss_df['chromEnd'], errors='coerce')
  tss_df = tss_df.dropna()
  tss_df['chromStart'] = tss_df['chromStart'].astype(int)
  tss_df['chromEnd'] = tss_df['chromEnd'].astype(int)
  tss_df = tss_df.sort_values(by=['#chrom', 'chromStart'])

  return tss_df


def get_upstream_mutations(chunk, tss_df, upstream):
  """
  sort and preprocess bed file containing transcription start sites
  :param chunk: 10 000 000 rows of somatic mutations in a pandas dataframe format
  :param tss_df: dataframe of transcription start site positions, sorted on chromosome position
  :param upstream: int representating how many bases to check upstream of TSS
  :return: a processed chunk of a pandas dataframe, containing only the somatic mutations within 1200 bases upstream of the TSS
  """
  upstream_chunk = pd.DataFrame()
  for chromosome in set(chunk['chromosome']):
    chr_chunk = chunk[chunk['chromosome'] == chromosome]
    # same chrosome format required
    chr_tss_df = tss_df[tss_df['#chrom'] == 'chr'+str(chromosome)]
    upstream_chr_chunk = pd.DataFrame()

    for tss_start in chr_tss_df['chromStart']:
      # only keep mutations within N amount of bases upstream of TSS
      tmp_chr_chunk = chr_chunk[chr_chunk['chromosome_start'].between(tss_start - upstream, tss_start, inclusive='left')]
      # skip concat attempt for increased speed
      if len(tmp_chr_chunk) == 0:
        continue
      upstream_chr_chunk = pd.concat([upstream_chr_chunk, tmp_chr_chunk])

    upstream_chunk = pd.concat([upstream_chunk, upstream_chr_chunk])

  return upstream_chunk


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-t", "--TSS", type=str, required=True, help="Path to BED of Transcription Start Sites of interest")
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Path to export the gene expression to")
  args = parser.parse_args()
  
  main(args)
