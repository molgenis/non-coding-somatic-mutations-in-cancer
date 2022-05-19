#!/usr/bin/env python3
"""
This script can be used to convert the .txt based eQTL data from eQTLGen to a valid BED format.
TODO: use same formatting as other scripts for meta data

./preprocess_eqtl_txt2bed.py -i ../../GREEN_DB/2019-12-11-cis-eQTLsFDR-Probe
Level-CohortInfoRemoved-BonferroniAdded.txt.gz -o ../../GREEN_DB/2022-05-13_cis-eqtl_above_p_0_9.bed -m -l -c
"""
__author__ = "Tijs van Lieshout"
# Imports
import argparse

from pybedtools import BedTool as bt
import pandas as pd
from tqdm import tqdm

def main(args):
  convert_txt2bed(args.inputPath, args.outputPath, args.outputMinimal, 
                  args.sortLexographically, args.controlSNPs)


def convert_txt2bed(input_path, output_path, output_minimal, sort_lexo, isControlSet):
  if input_path.endswith(".gz"):
    if isControlSet:
      chunk_size = 10 ** 6
      input_df = pd.DataFrame()
      for chunk in tqdm(pd.read_csv(input_path, sep="\t", compression='gzip', chunksize=chunk_size)):
        processed_chunk = get_control_snps(chunk)
        input_df = pd.concat([input_df, processed_chunk])
    else:
      input_df = pd.read_csv(input_path, sep="\t", compression='gzip')
  else:
    input_df = pd.read_csv(input_path, sep="\t")

  if output_minimal: # create the bare minimal BED format
    minimal_df = pd.DataFrame()
    minimal_df['chrom'] = input_df['SNPChr']
    minimal_df['chromStart'] = input_df['SNPPos']
    minimal_df['chromEnd'] = input_df['SNPPos']
    minimal_df['name'] = input_df['SNP']

    # sorting is required for many bedtools steps
    if sort_lexo: # default for BedTools 
      minimal_df['chrom'] = 'chr'+minimal_df['chrom'].astype(str)
      minimal_df = minimal_df.sort_values(by=['chrom', 'chromStart'])
    else: # default for this script if user has not supplied lexographically
      minimal_df = minimal_df.sort_values(by=['chrom', 'chromStart']) # sorting is required for many bedtools steps
      minimal_df['chrom'] = 'chr'+minimal_df['chrom'].astype(str)

    minimal_df = minimal_df.drop_duplicates()
    input_bt = bt.from_dataframe(minimal_df.rename(columns={'chrom': '#chrom'}), header=True)

  else:
    input_bt = bt.from_dataframe(input_df)
    print("immediate, and thus rough, conversion")

  input_bt.saveas(output_path)


def get_control_snps(chunk):
  chunk = chunk[chunk['Pvalue'] == 1]
  return chunk
  

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for txt to convert") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Output path for resulting BED file")
  parser.add_argument("-m", "--outputMinimal", action="store_true", required=False, help="Only keep chrom, start, stop and name")
  parser.add_argument("-l", "--sortLexographically", action="store_true", required=False, 
    help="Sort BED file lexographically, such that chr10 comes before chr2. This is the default behaviour for BedTools")
  parser.add_argument("-c", "--controlSNPs", action="store_true", required=False, 
    help="Specify if input are control SNPs and you want to only keep p-values near 1")
  args = parser.parse_args()
  
  main(args)

