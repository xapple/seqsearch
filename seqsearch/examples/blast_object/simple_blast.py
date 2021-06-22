#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

A script to demonstrate the usage of the `BLASTquery` class.
"""

# Built-in modules #
import inspect, os

# Internal modules #
from seqsearch.databases.ncbi_16s import ncbi_16s
from seqsearch.search.blast       import BLASTquery

# First party modules #
from fasta import FASTA

# Get current directory #
file_name = inspect.getframeinfo(inspect.currentframe()).filename
this_dir  = os.path.dirname(os.path.abspath(file_name)) + '/'

###############################################################################
if __name__ == "__main__":

    # Main input #
    seqs = FASTA(this_dir + 'seqs.fasta')

    # The database to search against #
    db = ncbi_16s.blast_db

    # Create search #
    query = BLASTquery(seqs, db)

    # Run #
    query.run()

