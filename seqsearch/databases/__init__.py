#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, fnmatch
from collections import Counter
import urllib.request

# First party modules #
from fasta import FASTA
from autopaths.auto_paths import AutoPaths
from autopaths.dir_path   import DirectoryPath
from autopaths.file_path  import FilePath
from plumbing.cache       import property_cached
from plumbing.common      import natural_sort
from plumbing.scraping    import download_from_url

# Third party modules #
from tqdm import tqdm

# Constants #
home = os.environ.get('HOME', '~') + '/'
base_directory = DirectoryPath(home + "databases/")

###############################################################################
class Database:
    """General database object to inherit from."""

    all_paths = """
    /raw/
    /unzipped/
    /blast_db/
    """

    def __init__(self, seq_type='nucl', base_dir=None):
        # The sequence type is either 'prot' or 'nucl' #
        self.seq_type = seq_type
        # The default base directory #
        if base_dir is None: base_dir = base_directory
        # Make base_dir object #
        self.base_dir = base_dir + self.short_name + '/'
        # Make autopaths object #
        self.autopaths = AutoPaths(self.base_dir, self.all_paths)

    def __repr__(self):
        # Get the name of this class #
        name = self.__class__.__name__
        # Return a user-friendly string #
        return '<%s at "%s">' % (name, self.base_dir)

    def __bool__(self):
        """
        Return True if the database was already downloaded and the
        results are stored on the filesystem. Return False otherwise.
        """
        return not self.autopaths.unzipped_dir.empty

    @property_cached
    def files_to_retrieve(self):
        """The files we want to download with their destinations."""
        return {f: FilePath(self.autopaths.raw_dir + f)
                for f in self.files}

    #---------------------------- GZIP compression ---------------------------#
    def ungzip(self):
        """Ungzip them."""
        # Check the extension #
        for f in tqdm(self.raw_files):
            if f.endswith('.gz') or f.endswith('.gzip'): continue
            destination = self.autopaths.unzipped_dir + f.prefix
            f.ungzip_to(destination)
        # Make them only readable #
        for f in self.autopaths.unzipped_dir:
            f.permissions.only_readable()

    def untargz(self):
        """Untargzip them."""
        # Check the extension #
        for f in tqdm(self.raw_files):
            if f.endswith('.tar.gz') or f.endswith('.tgz'): continue
            f.untargz_to(self.autopaths.unzipped_dir)
        # Make them only readable #
        for f in self.autopaths.unzipped_dir:
            f.permissions.only_readable()

    @property
    def raw_files(self):
        """The files we have downloaded."""
        return self.autopaths.raw_dir.contents

    @property
    def unzipped_files(self):
        """The files we have downloaded."""
        return self.autopaths.unzipped_dir.contents

    #------------------------ Only for FASTA databases -----------------------#
    @property
    def fasta_files(self):
        """The files we have downloaded."""
        return map(FASTA, self.raw_files)

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

    #-------------------- Only for BWA indexed databases ---------------------#
    @property
    def bwa_index(self):
        # Get the first file #
        first = next(self.unzipped_files)
        # Get only the prefix without the extension #
        return first.prefix_path

###############################################################################
class DatabaseHTTP(Database):
    """A database that is stored on an HTTP server."""

    @property_cached
    def files_to_retrieve(self):
        """The files we want to download with their destinations."""
        return {f: FilePath(self.autopaths.raw_dir + f)
                for f in self.files}

    @property
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        # Function to get the size of a file #
        def get_size_http(url):
            response = urllib.request.urlopen(url)
            return int(response.getheader("Content-Length"))
        # Check each file #
        return {source: dest
                for source, dest in self.files_to_retrieve.items()
                if dest.count_bytes == 0 or
                dest.count_bytes != get_size_http(self.base_url + source)}

    def download(self):
        """Retrieve all files from the website."""
        for source, dest in self.files_remaining.items():
            # Get the full URL #
            url = self.base_url + source
            # Similar to wget #
            download_from_url(url, dest,
                              stream   = True,
                              progress = True,
                              desc     = source,
                              cleanup  = True,
                              )
            # Make it readable only #
            dest.permissions.only_readable()

###############################################################################
class DatabaseFTP(Database):
    """A database that is stored on an FTP server."""

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
        # In the case we specify a pattern #
        if hasattr(self, "pattern"):
            files = self.ftp.listdir(self.ftp.curdir)
            files.sort(key=natural_sort)
            return {f: FilePath(self.autopaths.raw_dir + f)
                    for f in files if fnmatch.fnmatch(f, self.pattern)}
        # In the case we specify a list of files #
        if hasattr(self, "files"):
            return {f: FilePath(self.autopaths.raw_dir + f)
                    for f in self.files}

    @property
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        return {source: dest for source, dest in self.files_to_retrieve.items()
                if dest.count_bytes != self.ftp.path.getsize(source)}

    def download(self):
        """Retrieve all files from the FTP site."""
        # Create the directory #
        self.base_dir.create_if_not_exists()
        # Loop over files #
        for source, dest in tqdm(self.files_remaining.items()):
            dest.remove()
            self.ftp.download(source, dest)
            dest.permissions.only_readable()