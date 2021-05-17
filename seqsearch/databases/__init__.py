#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, fnmatch
from collections import OrderedDict, Counter

# First party modules #
from fasta import FASTA
from autopaths.auto_paths import AutoPaths
from autopaths.dir_path   import DirectoryPath
from autopaths.file_path  import FilePath
from plumbing.cache       import property_cached
from plumbing.common      import natural_sort

# Third party modules #
from tqdm import tqdm

# Constants #
home = os.environ.get('HOME', '~') + '/'
base_directory = home + "databases/"

###############################################################################
class Database:
    """General database object to inherit from."""

    all_paths = """
    /raw/
    /unzipped/
    /blast_db/
    """

    def __init__(self, seq_type=None, base_dir=None):
        # The sequence type is either 'prot' or 'nucl' #
        self.seq_type = seq_type
        # The default base directory #
        if base_dir is None:
            base_dir = os.environ.get('HOME', '/') + '/'
        # Make base_dir object #
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.base_dir = DirectoryPath(self.base_dir)
        # Make autopaths object #
        self.autopaths = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def ftp(self):
        """If the data is to be obtained by FTP, here is the ftputil object."""
        from ftputil import FTPHost
        ftp = FTPHost(self.ftp_url, "anonymous")
        ftp.chdir(self.ftp_dir)
        return ftp

    @property_cached
    def files_to_retrieve(self):
        """The files we want to download with their destinations."""
        if hasattr(self, "pattern"):
            files = self.ftp.listdir(self.ftp.curdir)
            files.sort(key=natural_sort)
            return OrderedDict((f, FilePath(self.autopaths.raw_dir+f)) for f in files
                               if fnmatch.fnmatch(f, self.pattern))
        if hasattr(self, "files"):
            return OrderedDict((f, FilePath(self.autopaths.raw_dir+f)) for f in self.files)

    @property
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        return OrderedDict((source,dest) for source, dest in self.files_to_retrieve.items()
                           if dest.count_bytes != self.ftp.path.getsize(source))

    def download(self):
        """Retrieve all files from the FTP site."""
        # Create the directory #
        self.base_dir.create_if_not_exists()
        # Loop over files #
        for source, dest in tqdm(self.files_remaining.items()):
            dest.remove()
            self.ftp.download(source, dest)
            dest.permissions.only_readable()

    @property
    def raw_files(self):
        """The files we have downloaded."""
        return map(FASTA, self.autopaths.raw_dir.contents)

    def ungzip(self):
        """Ungzip them."""
        # Gzip #
        for f in tqdm(self.raw_files):
            destination = self.autopaths.unzipped_dir+f.prefix
            f.ungzip_to(destination)
            destination.permissions.only_readable()

    def untargz(self):
        """Untargzip them."""
        # Gzip #
        for f in tqdm(self.raw_files): f.untargz_to(self.autopaths.unzipped_dir)
        for f in self.autopaths.unzipped_dir: f.permissions.only_readable()

    @property
    def sequences(self):
        """All the sequences from all the raw files."""
        for fasta in self.raw_files:
            for seq in fasta: yield seq

    #------------------ Only for preformatted BLAST databases ----------------#
    @property_cached
    def blast_db(self):
        """A BLASTable version of the sequences."""
        # Import #
        from seqsearch.search.blast import BLASTdb
        # Create object #
        db = BLASTdb(self.autopaths.unzipped_dir + self.db_name,
                     self.seq_type)
        # Return #
        return db

    #--------------------- Only for taxonomic databases ----------------------#
    @property_cached
    def tax_depth_freq(self):
        def depths():
            with open(self.taxonomy, 'r') as handle:
                for line in handle:
                    line = line.strip('\n')
                    otu_name, species = line.split('\t')
                    yield len(species.split(';'))
        return Counter(depths())
