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
from fasta import FASTA
from autopaths.auto_paths import AutoPaths

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class Silva(Database):
    """
    SILVA provides comprehensive, quality checked and regularly updated
    datasets of aligned small (16S/18S, SSU) and large subunit (23S/28S, LSU)
    ribosomal RNA (rRNA) sequences for all three domains of life
    (Bacteria, Archaea and Eukarya).
    SILVA are the official databases of the software package ARB.

    https://www.arb-silva.de

    To install:
        from seqsearch.databases.silva import silva
        silva.download()
        silva.unzip()

    It will put it in ~/databases/silva_xxx/
    """

    view_url   = "https://www.arb-silva.de/no_cache/download/archive/"
    base_url   = "https://www.arb-silva.de/fileadmin/silva_databases/"
    short_name = "silva"

    all_paths = """
    /test.txt
    """

    def __init__(self, version, seq_type, base_dir=None):
        # Attributes #
        self.version    = version
        self.seq_type   = seq_type
        self.short_name = self.short_name + "_" + self.version
        # Base directory #
        if base_dir is None: base_dir = home
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # URL #
        self.url  = "release_%s/Exports/"  % self.version
        # The database #
        self.nr99_name = "SILVA_%s_SSURef_Nr99_tax_silva.fasta.gz" % self.version
        self.nr99_dest = FASTA(self.base_dir + self.nr99_name)
        self.nr99      = FASTA(self.base_dir + self.nr99_name[:-3])
        # The alignment #
        self.aligned_name = "SILVA_%s_SSURef_Nr99_tax_silva_full_align_trunc.fasta.gz" % self.version
        self.aligned_dest = FASTA(self.base_dir + self.aligned_name)
        self.aligned      = FASTA(self.base_dir + self.aligned_name[:-3])

    def download(self):
        self.nr99_dest.directory.create(safe=True)
        self.nr99_dest.remove(safe=True)
        self.aligned_dest.remove(safe=True)
        import wget
        print("\nDownloading", self.base_url + self.url + self.nr99_name)
        wget.download(self.base_url + self.url + self.nr99_name,    out=self.nr99_dest.path)
        print("\nDownloading", self.base_url + self.url + self.aligned_name)
        wget.download(self.base_url + self.url + self.aligned_name, out=self.aligned_dest.path)

    def unzip(self):
        self.nr99_dest.ungzip_to(self.nr99)
        self.nr99.permissions.only_readable()
        self.aligned_dest.ungzip_to(self.aligned)
        self.aligned.permissions.only_readable()

###############################################################################
silva = Silva("128", "nucl")