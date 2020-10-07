#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os

# First party modules #
from seqsearch.databases import Database

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class NCBI16S(Database):
    """
    The NCBI provides a database of 16S ribosomal sequences.
    The exact title is "Bacteria and Archaea: 16S ribosomal RNA project".
    More information is available at:

    https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/

    To install:

        >>> from seqsearch.databases.ncbi_16s import ncbi_16s
        >>> ncbi_16s.download()
        >>> ncbi_16s.untargz()

    It will place the resulting files in "~/databases/ncbi_16s/".
    """

    tag        = "ncbi_16s"
    short_name = "ncbi_16s"
    long_name  = 'The NCBI 16S RNA database'

    base_url   = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/blast/db/"
    files      = ['16S_ribosomal_RNA.tar.gz']
    db_name    = '16S_ribosomal_RNA'

###############################################################################
ncbi_16s = NCBI16S('nucl')