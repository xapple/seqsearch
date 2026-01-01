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
from seqsearch.databases import DatabaseFTP

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class HumanGenome(DatabaseFTP):
    """
    The NCBI provides the latest version of the human genome, as well
    as preformatted files here:

    https://www.ncbi.nlm.nih.gov/genome/guide/human/

    To install:

        >>> from seqsearch.databases.human import hg38
        >>> hg38.download()
        >>> hg38.untargz()
        >>> hg38.autopaths.raw_dir.remove()

    It will place the resulting files in "~/databases/human/".
    """

    tag        = "hg38"
    short_name = "grc_h38"
    long_name  = 'Human Genome v38 at NCBI'

    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/" \
                 "seqs_for_alignment_pipelines.ucsc_ids/"

    files = ['GCA_000001405.15_GRCh38_full_analysis_set.fna.bwa_index.tar.gz']

###############################################################################
# Create a singleton #
hg38 = HumanGenome()