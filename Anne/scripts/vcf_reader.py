#!/usr/bin/env python3

# Imports
import pandas as pd
import io


def read_vcf(path):
    """
    Read vcf files.
    Comes from: https://gist.github.com/dceoy/99d976a2c01e7f0ba1c813778f9db744
        Author: Daichi Narushima (成島 大智)
        GitName: dceoy
        File: read_vcf.py
    :param path: Path to the file
    :return:
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
