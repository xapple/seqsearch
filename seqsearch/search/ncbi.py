#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import time

# Internal modules #

# First party modules #
from tqdm import tqdm

# Third party modules #
from Bio import Entrez
from Bio.Entrez.Parser import CorruptedXMLError

# Constants #
Entrez.email = "I don't know who will be running this script"

###############################################################################
class UtilsNCBI(object):
    """
    An object that takes care of extracting NCBI information via the eutils
    http://www.ncbi.nlm.nih.gov/books/NBK25499/
    """

    def gi_num_to_tax(self, id_num):
        """
        How to convert a single GI identification number to taxonomy info
        Can also accept a list of GI numbers in the parameter `id_num`
        """
        if isinstance(id_num, list):
            gb_entries = Entrez.efetch(db="nuccore", id=id_num, rettype="fasta", retmode="xml")
            gb_records = Entrez.read(gb_entries)
            return gb_records
        else:
            gb_entry = Entrez.efetch(db="nuccore", id=id_num, rettype="fasta", retmode="xml")
            gb_records = Entrez.read(gb_entry)
            tax_num = gb_records[0]['TSeq_taxid']
            handle = Entrez.efetch(db="Taxonomy", id=tax_num, retmode="xml")
            records = Entrez.read(handle)
            return records[0]['Lineage']

    #-------------------------------------------------------------------------#
    def gis_to_records(self, gis, progress=True):
        """Download information from NCBI in batch mode.
        Return a dictionary with gi numbers as keys and records
        as values."""
        # Should we display progress ? #
        progress = tqdm if progress else lambda x:x
        # Do it by chunks #
        gis = list(gis)
        at_a_time = 400
        result = {}
        # Main loop #
        for i in progress(range(0, len(gis), at_a_time)):
            chunk = gis[i:i+at_a_time]
            chunk = map(str, chunk)
            records = self.chunk_to_records(chunk)
            result.update(dict(zip(chunk, records)))
        # Return #
        return result

    def chunk_to_records(self, chunk, validate=False):
        """
        Download from NCBI until it works. Will restart until reaching the python
        recursion limit. We don't want to get banned from NCBI so we have a little
        pause at every function call.
        """
        from urllib2 import HTTPError
        time.sleep(0.5)
        try:
            response = Entrez.efetch(db="nuccore", id=chunk, retmode="xml")
            records  = list(Entrez.parse(response, validate=validate))
            return records
        except (HTTPError, CorruptedXMLError):
            print("\nFailed downloading %i records, trying again\n" % len(chunk))
            return self.chunk_to_records(chunk)

    #-------------------------------------------------------------------------#
    def record_to_taxonomy(self, record): return record['GBSeq_taxonomy']

    def record_to_source(self, record):
        qualifiers = record['GBSeq_feature-table'][0]['GBFeature_quals']
        for qualifier in qualifiers:
            if qualifier['GBQualifier_name'] == 'isolation_source':
                return qualifier['GBQualifier_value']
