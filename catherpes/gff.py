#!/usr/bin/env python3

"""The catherpes gff.py module provides a class and methods for
parsing GFF3 formatted text files.  GFF3 is a file format 

Example:
    Give a few examples of using the module::

        $ python gff.py

Just because you are curious: Catherpes is the genus name of the
Canyon Wren found throughout the Western US and Mexico.  It is
indelibly linked in my mind with floating through the river canyons of
the American Southwest.  This library is named as a tribute to this
little bird, because, if I hear the distinctive and beautiful call of
the Canyon Wren, then I'm almost certainly in a very happy place!

Todo:
    * Write the damn code...
    * You have to also use ``sphinx.ext.todo`` extension

"""

__author__ = "Barry Moore"
__version__ = "0.1.0"
__license__ = "GNU GPL"

import argparse
import gzip
import numpy as np
import pandas as pd

def main(args):
    """ Main entry point of the app """
    print("catherpes/GFF")
    print(args)

    gff = GFF(file=args.file)
    headers = gff.headers[0:20]
    data = gff.data[0:200]

    for header in headers:
        print(header)

    for record in data:
        keys = [ 'seqid', 'attributes']
        print(record)
        # print({k:record[k] for k in keys if k in record})

class GFF(object):
    """Catherpes GFF is a Python class with methods for parsing and
    manipulating GFF3 data.
    """

    def __init__(self, file=None, format='dict'):
        """Args:
            file (str)  : The path/name of the GFF3 file to parse.

            format (str): The format to return the data as when the
                          data attribute is called.  Valid values are:

                          dict: returns the data as a
                                list-of-dictionaries and is the
                                default.

                          list: returns the data as a list-of-lists

                          object: returns the data as Catherpes
                                  Feature objects. Not implimented
                                  yet.

                          df: returns the data as a pandas
                                  dataframe. Not implimented yet.

                          xarray: returns the data as an Xarray.  Not
                                  implimented yet.

        """

        # Define attributes
        self.file = file
        self.format = format
        self.headers = []
        self.data = []

        # Parse file
        self._parse(file=file)

        if self.format == 'pandas':
            self.data = pd.DataFrame(self.data)

    def _parse(self, file=None):
        """
        Parse a GFF3 file.

        Args:
            file: The path/name of the GFF3 file to parse.

        Returns:
            A catherpes/GFF3 object.
        """

        keys = ('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes') 
        
        with gzip.open(file, 'rt') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('#'):
                    self.headers.append(line)
                else:
                    values = line.split('\t')
                    record = dict(zip(keys, values))
                    record['attributes'] = self._parse_attributes(record['attributes'])
                    self._promote_attributes(record)
                    self.data.append(record)

    def _parse_attributes(self, attrb_text):
        """Parse attribures in a GFF3 record.

        Args:
            attrbs: A string of text in GFF3 attributes column format.
                    This format is basically
                    'key1=value1;key2=value2;key3,key4'.

        Returns:
            A dictionary of attributes.

        """
        
        attrbs = {}
        for pair in attrb_text.split(';'):
            (key, value) = pair.split('=')
            value = True if value is None else value
            attrbs[key] = value

        return attrbs

    def _promote_attributes(self, record):
        """Promote a key subset of GFF3 attributes to primary level
        keys/columns in the parsed data structure.

        Args: None

        Returns: No return value.
        """

        prm_attrbs = ('ID', 'Name', 'Alias', 'Parent')

        for attr in prm_attrbs:
            if attr in record['attributes']:
                record[attr] = record['attributes'][attr]
            else:
                record[attr] = None

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("file", help="Required path/name of a GFF3 file")

    # Optional argument flag which defaults to False
    parser.add_argument("-f", "--flag", action="store_true", default=False)

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-n", "--name", action="store", dest="name")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
