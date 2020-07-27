#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from seqsearch.databases import Database

###############################################################################
class RefSeqBacteriaProtNR(Database):
    """
    The RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast database.
    """

    short_name = "refseq_bact_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/bacteria/"
    pattern    = 'bacteria.nonredundant_protein.*.protein.faa.gz'

###############################################################################
class RefSeqArchaeaProtNR(Database):
    """
    The RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast database.
    """

    short_name = "refseq_arch_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/archaea/"
    pattern    = 'archaea.nonredundant_protein.*.protein.faa.gz'

###############################################################################
refseq_bact_prot_nr = RefSeqBacteriaProtNR('prot')
refseq_arch_prot_nr = RefSeqArchaeaProtNR('prot')