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
from seqsearch.databases import base_directory

# First party modules #
from fasta import FASTA
from autopaths.auto_paths import AutoPaths
from autopaths.file_path import FilePath

###############################################################################
class String(object):
    """
    The STRING database. See:
    http://string.embl.de/newstring_cgi/show_download_page.pl
    """

    base_url = "http://string.embl.de/newstring_download/"
    short_name = "string"

    all_paths = """
    /raw/all_proteins.fasta.gz
    /raw/cog_mappings.tsv.gz
    /unzipped/all_proteins.fasta
    /unzipped/cog_mappings.tsv
    /blast_db/all_proteins.fasta
    /blast_db/all_proteins.fasta.00.pin
    /blast_db/logfile.txt
    /blast_db/out.txt
    """

    def __init__(self, seq_type='prot'):
        self.seq_type = seq_type
        self.base_dir = base_directory + self.short_name
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property
    def files_to_retrieve(self):
        """The files we want to download with their destinations."""
        result = OrderedDict()
        result[self.base_url + "protein.sequences.v9.1.fa.gz"] = FilePath(self.p.raw_proteins)
        result[self.base_url + "COG.mappings.v9.1.txt.gz"]     = FilePath(self.p.raw_mappings)
        return result

    @property
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        get_size_http = lambda url: urllib.urlopen(url).info().getheaders("Content-Length")[0]
        return OrderedDict((source, dest) for source, dest in self.files_to_retrieve.items()
                           if dest.count_bytes != get_size_http(source))

    def download(self):
        """Retrieve all files from the website"""
        for source, dest in self.files_remaining.items():
            dest.remove()
            urllib.urlretrieve(source, dest)
            dest.permissions.only_readable()

    @property
    def raw_files(self):
        """The files we have downloaded."""
        return map(FilePath, self.p.raw_dir.contents)

    def unzip(self):
        """Unzip them"""
        for f in self.raw_files: f.ungzip_to(self.p.unzipped_dir + f.prefix)

    @property
    def all_proteins(self):
        """The main fasta file."""
        return FASTA(self.p.unzipped_proteins)

    @property
    def mappings(self):
        """The cog mappings."""
        return FilePath(self.p.unzipped_mappings)

    @property
    def blast_db(self):
        """A BLASTable version of the sequences."""
        if not self.p.blast_fasta.exists:
            self.p.unzipped_proteins.link_to(self.p.blast_fasta, safe=True)
        from seqsearch.search.blast import BLASTdb
        blast_db = BLASTdb(self.p.blast_fasta, 'prot')
        if not self.p.pin.exists:
            blast_db.makeblastdb(logfile=self.p.logfile, out=self.p.out)
        return blast_db

###############################################################################
string = String()