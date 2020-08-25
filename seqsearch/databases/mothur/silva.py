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
class SilvaMothur(Database):
    """
    This is the SILVA database, in its specific version from mothur. Seen at:

    https://www.mothur.org/wiki/Silva_reference_files

    To install:

        >>> from seqsearch.databases.mothur.silva import silva_mothur
        >>> silva_mothur.download()
        >>> silva_mothur.unzip()

    It will place the results in `~/databases/silva_mothur/`.
    """

    tag        = "silva"
    short_name = "silva_mothur"
    long_name  = "The Silva v138 database (mothur version)"
    version    = "138"
    base_url   = "https://mothur.s3.us-east-2.amazonaws.com/wiki/"

    all_paths = """
                /silva.tgz
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
        self.url = self.base_url + "silva.nr_v%s.tgz" % self.version
        # Location of zip file locally #
        self.dest = self.autopaths.tgz
        # The results after download #
        self.alignment = self.base_dir + "silva.nr_v%s.align"
        self.taxonomy  = self.base_dir + "silva.nr_v%s.tax"
        # Make them FilePaths objects #
        self.alignment = FilePath(self.alignment % self.version)
        self.taxonomy  = FilePath(self.taxonomy  % self.version)
        # The part that mothur will use for naming files #
        self.nickname = "nr_v%s" % self.version

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
silva_mothur = SilvaMothur()