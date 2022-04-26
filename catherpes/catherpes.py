#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Barry Moore"
__version__ = "0.1.0"
__license__ = "GNU GPL"

import sys
import argparse
import click

@click.command()
def main(args=None):
    """Console script for catherpes."""
    click.echo("Replace this message by putting your code into "
               "catherpes.cli.main")
    click.echo("See click documentation at https://click.palletsprojects.com/")
    return 0

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
