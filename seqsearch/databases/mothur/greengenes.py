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
from autopaths.dir_path   import DirectoryPath

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class GreengenesMothur(Database):
    """
    This is the Greengenes database, in its specific version from mothur.
    Seen at:

    https://mothur.org/wiki/greengenes-formatted_databases/

    To install:

        >>> from seqsearch.databases.mothur.greengenes import gg_mothur
        >>> gg_mothur.download()
        >>> gg_mothur.unzip()

    It will place the results in `~/databases/gg_mothur/`.

    This database is from 2013.
    """

    nickname   = "gg"
    tag        = "greengenes"
    short_name = "gg_mothur"
    long_name  = "The Greengenes v13_8_99 database (mothur version)"
    version    = "13_8_99"
    base_url   = "https://mothur.s3.us-east-2.amazonaws.com/wiki/"

    all_paths = """
                /gg_alignment.tgz
                /gg_taxonomy.tgz
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
        self.base_dir  = DirectoryPath(data_dir + self.short_name + '/')
        self.autopaths = AutoPaths(self.base_dir, self.all_paths)
        # Location of zip file remotely #
        self.ref_url = self.base_url + "gg_13_8_99.refalign.tgz"
        self.tax_url = self.base_url + "gg_13_8_99.taxonomy.tgz"
        # Location of zip file locally #
        self.ref_dest = self.autopaths.alignment
        self.tax_dest = self.autopaths.taxonomy
        # The results after download #
        self.alignment = self.base_dir + "gg_13_8_99.refalign"
        self.taxonomy  = self.base_dir + "gg_13_8_99.gg.tax"
        # Make them FilePaths objects #
        self.alignment = FilePath(self.alignment)
        self.taxonomy  = FilePath(self.taxonomy)

    def download(self):
        # Make sure the directory exists #
        self.base_dir.create(safe=True)
        # Remove previous downloads #
        self.ref_dest.remove()
        self.tax_dest.remove()
        # Message #
        print("\n Downloading '%s'" % self.ref_url)
        # Download #
        import wget
        wget.download(self.ref_url, out=self.ref_dest.path)
        # Message #
        print("\n Downloading '%s'" % self.tax_url)
        # Download #
        wget.download(self.tax_url, out=self.tax_dest.path)

    def unzip(self):
        # Message #
        print("\n Extracting archive '%s'" % self.ref_dest)
        # Uncompress #
        archive = tarfile.open(self.ref_dest, 'r:gz')
        archive.extractall(self.base_dir)
        # Message #
        print("\n Extracting archive '%s'" % self.tax_dest)
        # Uncompress #
        archive = tarfile.open(self.tax_dest, 'r:gz')
        archive.extractall(self.base_dir)

    def __bool__(self):
        """
        Return True if the silva database was already downloaded and the
        results are stored on the filesystem. Return False otherwise.
        """
        return self.taxonomy.exists and self.alignment.exists

###############################################################################
gg_mothur = GreengenesMothur()