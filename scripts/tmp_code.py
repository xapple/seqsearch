#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio

Development script to test some of the methods in the `seqsearch` package
and try out different things. This script can safely be ignored
and is meant simply as a sandbox.

Typically you would run this file from a command line like this:

     ipython3 -i -- ~/deploy/seqsearch/scripts/tmp_code.py
"""

# Built-in modules #

# Internal modules #
from seqsearch.databases.ncbi_16s import ncbi_16s

# Third party modules #

# Constants #

###############################################################################
# Create project #
ncbi_16s.download()
ncbi_16s.untargz()