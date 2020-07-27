#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import urllib
from collections import OrderedDict

# Internal modules #
from seqsearch.search.blast import BLASTdb
from seqsearch.databases import base_directory

# First party modules #
from fasta import FASTA
from autopaths.auto_paths import AutoPaths
from autopaths.file_path import FilePath
from plumbing.csv_tables import TSVTable

###############################################################################
class NonRedundant(object):
    """
    The NR database from NCBI.
    NR contains non-redundant sequences from GenBank translations
    (i.e. GenPept) together with sequences from other databanks
    (Refseq, PDB, SwissProt, PIR and PRF).
    """

