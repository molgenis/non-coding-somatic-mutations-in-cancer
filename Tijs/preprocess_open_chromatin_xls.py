#!/usr/bin/env python3
"""
Script used for processing a excel file of a bed file per sheet on chromatin states

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
./preprocess_open_chromatin_xls.py -i ../../chromatin_states/2022-05-12_1256271tables1.xls
"""

# Metadata
__title__ = "Process a excel file of a bed file per sheet on chromatin states"
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-12"
__updated__ = "2022-05-12"
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


def main(args):
  sheets = pd.read_excel(args.inputPath, sheet_name=None)

  for name, df in sheets.items():
    df = df.iloc[:, :4]
    df = df.rename(columns={df.columns[0]: 'name', df.columns[1]: '#chr'})
    df = df[[df.columns[1], df.columns[2], df.columns[3], df.columns[0]]] 
    
    bt.from_dataframe(df, header=True).saveas(f"{name}.bed")


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Input path for .xls to convert") 
  args = parser.parse_args()
  
  main(args)
