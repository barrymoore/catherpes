#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Created on Fri Jun 24 17:39:58 UTC 2022

Synopsis:

viq_batch_report --base viq_batch_report --manifest viq_output_manifest.

Description:

A script to parse multiple vIQ output files (typically in the dozens
or more because summary stats are provided).

The default behavior is to produce several summary files that
aggregate the data for various sections of a vIQ output

Args:

  --manifest, -m

    A two-column, tab-delimited file that provides an ID and vIQ file
    name for vIQ output to be included in the report(s).  This
    argument is required.  The two columns are:

        1) ID: A unique identifier for the vIQ report
        2) File name: The file name (with absolute or relative path) to 

  --base, -b [viq_batch_report]

    Base name for output files.

"""

import sys
import argparse
import json
import pandas as pd

sys.path.append('../')
from catherpes.viq import VIQ

def main():
    parser = argparse.ArgumentParser(
        description='Program description',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--manifest',
                        '-m',
                        help='A manifest file of ID\tFilename pairs for each vIQ output')
    parser.add_argument('--base',
                        '-b',
                        default='viq_batch_report',
                        help='The base name (prefix) to use for output files')
    parser.add_argument('--genes_only',
                        '-g',
                        default=False,
                        action="store_true",
                        help='Print only the gene records as tsv.')
    parser.add_argument('--meta_only',
                        '-e',
                        default=False,
                        action="store_true",
                        help='Print only the metadata as tsv.')
    args = parser.parse_args()

    fh = open(args.manifest, "r")

    samples = []
    for line in fh:
        line = line.strip()
        (sample_id, filename) = line.split("\t")
        samples.append((sample_id, filename))

    fh.close()
        
    if args.genes_only:
        print_genes(samples)
    elif args.meta_only:
        print_meta(samples)

def genes_df(samples):
    data_frames = []
    for (sample_id, filename) in samples:
        viq = VIQ(filename)
        df = viq.genes_as_df()
        df.insert(0, 'sample_id', sample_id)
        data_frames.append(df)
        
    return pd.concat(data_frames)

def print_genes(samples, filename='-'):
    df = genes_df(samples)
    if filename == '-':
        print(df.to_csv(sep='\t', index=False))
    else:
        df.to_csv(filename, sep='\t', index=False)

def meta_df(samples):
    data_frames = []
    for (sample_id, filename) in samples:
        viq = VIQ(filename)
        df = viq.meta_as_df()
        df.insert(0, 'sample_id', sample_id)
        data_frames.append(df)

    return pd.concat(data_frames)

def print_meta(samples, filename='-'):
    df = meta_df(samples)
    if filename == '-':
        print(df.to_csv(sep='\t', index=False))
    else:
        df.to_csv(filename, sep='\t', index=False)

if __name__ == "__main__":
    main()
