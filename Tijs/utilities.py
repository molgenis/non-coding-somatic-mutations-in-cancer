#!/usr/bin/env python3
"""
Utility functions for other scripts

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
NA
"""

# Metadata
__title__ = "Utility functions for other scripts" 
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-04"
__updated__ = "2022-05-04"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 1.0
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import yaml

import numpy as np
import pandas as pd
from pybedtools import BedTool as bt


def get_config():
    '''
    Open a yaml formatted config file in a safe way.
    '''
    with open('config.yaml', 'r') as stream:
        config = yaml.safe_load(stream)
    return config


def preprocess_bed_file(path, sort=False):
  bed_df = pd.read_csv(path, sep="\t")
  bed_df[bed_df.columns[0]] = np.where(bed_df[bed_df.columns[0]].str.contains("chr"), 
                                       bed_df[bed_df.columns[0]], 
                                       'chr' + bed_df[bed_df.columns[0]])
  bed = bt.from_dataframe(bed_df)

  if sort:
    bed = bed.sort()

  return bed

