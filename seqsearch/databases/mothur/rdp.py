#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, tarfile

# First party modules #
from seqsearch.databases  import Database
from autopaths.auto_paths import AutoPaths
from autopaths.file_path  import FilePath

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class RdpMothur(Database):
    """
    This is the RDP database, in its specific version from mothur. Seen at:

    https://mothur.org/wiki/rdp_reference_files/

    To install:

        >>> from seqsearch.databases.mothur.rdp import rdp_mothur
        >>> rdp_mothur.download()
        >>> rdp_mothur.unzip()

    It will place the results in `~/databases/rdp_mothur/`.

    This database is from March 2016.
    """

    nickname   = "rdp"
    tag        = "rdp"
    short_name = "rdp_mothur"
    long_name  = "The RDP v16 database (mothur version)"
    version    = "16"
    base_name  = "trainset16_022016.rdp"
    base_url   = "https://mothur.s3.us-east-2.amazonaws.com/wiki/"

    all_paths = """
                /rdp.tgz
                """

    @property
    def rank_names(self):
        """
        The names of the taxonomic rank at each level.
        There are a total of 7 ranks.
        """
        return ['Domain',    # 0
                'Phylum',    # 1
                'Class',     # 2
                'Order',     # 3
                'Family',    # 4
                'Genus',     # 5
                'Species']   # 6

    def __init__(self, data_dir=None):
        # The directory that contains all databases #
        if data_dir is None: data_dir = home + 'databases/'
        # Base directory for paths #
        self.base_dir  = data_dir + self.short_name + '/'
        self.autopaths = AutoPaths(self.base_dir, self.all_paths)
        # Location of zip file remotely #
        self.url = self.base_url + self.base_name + ".tgz"
        # Location of zip file locally #
        self.dest = self.autopaths.tgz
        # The results after download #
        prefix = self.base_dir + self.base_name + '/' + self.base_name
        self.alignment = FilePath(prefix + ".fasta")
        self.taxonomy  = FilePath(prefix + ".tax")

    def download(self):
        # Make sure the directory exists #
        self.dest.directory.create(safe=True)
        # Remove previous download #
        self.dest.remove()
        # Message #
        print("\n Downloading '%s'" % self.url)
        # Download #
        import wget
        return wget.download(self.url, out=self.dest.path)

    def unzip(self):
        # Message #
        print("\n Extracting archive '%s'" % self.dest)
        # Uncompress #
        archive = tarfile.open(self.dest, 'r:gz')
        archive.extractall(self.base_dir)

    def __bool__(self):
        """
        Return True if the silva database was already downloaded and the
        results are stored on the filesystem. Return False otherwise.
        """
        return self.taxonomy.exists

###############################################################################
rdp_mothur = RdpMothur()