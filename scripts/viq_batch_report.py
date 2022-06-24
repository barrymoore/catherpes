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

def main():
    import sys
    import argparse
    import json

    sys.path.append('../')
    from catherpes.viq import VIQ
    
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
    args = parser.parse_args()

    fh = open(args.manifest, "r")

    samples = []
    for line in fh:
        line = line.strip()
        (sample_id, filename) = line.split("\t")
        samples.append((sample_id, filename))
    fh.close()

    for (sample_id, filename) in samples:
        viq_data = VIQ(filename)
        
        for record in viq_data.gene_records.values():
            print("\t".join([sample_id, record.gene, record.viqscr, json.dumps(record.var_qual)]))        

def my_function():
    # Do something
    pass

if __name__ == "__main__":
    main()

