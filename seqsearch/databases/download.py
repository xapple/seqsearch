#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #

# Third party modules #

###############################################################################
def acc_to_fasta(accessions):
    """
    Pass a list of accessions IDs as argument and a string representing
    a FASTA is returned.
    """
    from Bio import Entrez
    Entrez.email = "I don't know who will be running this script"
    entries = Entrez.efetch(db      = "nuccore",
                            id      = accessions,
                            rettype = "fasta",
                            retmode = "xml")
    records = Entrez.read(entries)
    return records