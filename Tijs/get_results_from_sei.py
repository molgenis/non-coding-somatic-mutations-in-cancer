#!/usr/bin/env python3
"""
get top 3 diffs from DeepSEA-sei run

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
./get_results_from_sei.py -i ../../deepsea-sei/sei-framework/RPS26_UT_unique/c7ce95e5-6f9d-4249-a9bd-a01e3bda8689_RPS26_summary_file_UT_unique_diffs.tsv -o 2022-07-18_sei_top3_diffs.tsv
"""

# Metadata
__title__ = "get top 3 diffs from DeepSEA-sei run" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-07-18"
__updated__ = "2022-07-18"
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

from tqdm import tqdm


def main(args):
  df = pd.DataFrame()
  # TODO: include top_sequence_class
  seq_class_df = open_seq_class_df(args.inputPath)

  df = get_top3_sei(args.inputPath, seq_class_df)
  df = df.sort_values(by=["score1"], ascending=False)
  df.to_csv(args.outputPath, sep="\t", index=False)


def open_seq_class_df(path):
  split_path = path.split("_")
  split_path[-1] = "sequence-class-scores.tsv"
  seq_path = "_".join(split_path)

  return pd.read_csv(seq_path, sep="\t")


def get_top3_sei(input_path, seq_class_df):
  df = pd.DataFrame()
  for chunk in tqdm(pd.read_csv(input_path, sep="\t", chunksize=1)):
    chunk = chunk.squeeze()

    seq_class_entry = seq_class_df[(seq_class_df["chrom"] == chunk["chrom"]) & (seq_class_df["pos"] == chunk["position"])]
    for column in seq_class_entry.drop(columns=["seqclass_max_absdiff"]):
      if seq_class_entry[column].values == seq_class_entry["seqclass_max_absdiff"].values:
        seq_class = column
        break

    scores = chunk[8:].astype(np.float32)
    top3 = scores.nlargest(3)
    entry = pd.DataFrame.from_dict({"chr":[chunk["chrom"]], "position":[chunk["position"]], "top1":[top3.index[0]], "top2":[top3.index[1]], "top3":[top3.index[2]], 
                                    "score1":[top3[0]], "score2":[top3[1]], "score3":[top3[2]], "top_sequence_class":[seq_class], "top_sequence_class_score":[seq_class_entry["seqclass_max_absdiff"].values[0]]})
    df = pd.concat([df, entry])

  return df


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  
  main(args)
