# Built-in modules #
import os, fnmatch
from collections import OrderedDict

# First party modules #
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.cache import property_cached
from plumbing.common import natural_sort
from fasta import FASTA

# Third party modules #
from ftputil import FTPHost
from tqdm import tqdm

# Constants #
home = os.environ['HOME'] + '/'
base_directory = home + "/databases/"

###############################################################################
class Database(object):
    """General database object to inherit from."""

    all_paths = """
    /raw/
    /blast_db/
    """

    def __init__(self, seq_type):
        # Attributes #
        self.seq_type = seq_type
        self.base_dir = base_directory + self.short_name
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def ftp(self):
        """If the data is to be obtained by FTP here is the ftputil object."""
        ftp = FTPHost(self.ftp_url, "anonymous")
        ftp.chdir(self.ftp_dir)
        return ftp

    @property_cached
    def files_to_retrive(self):
        """The files we want to download with their destinations."""
        files = self.ftp.listdir(self.ftp.curdir)
        files.sort(key=natural_sort)
        return OrderedDict((f, FilePath(self.p.raw_dir+f)) for f in files
                            if fnmatch.fnmatch(f, self.pattern))

    @property
    def files_remaining(self):
        """The files we haven't downloaded yet based on size checks."""
        return OrderedDict((source,dest) for source,dest in self.files_to_retrive.items()
                           if dest.count_bytes != self.ftp.path.getsize(source))

    def download(self):
        """Retrieve all files from the FTP site"""
        for source,dest in tqdm(self.files_remaining.items()):
            dest.remove()
            self.ftp.download(source, dest)
            dest.permissions.only_readable()

    @property
    def raw_files(self):
        """The files we have downloaded."""
        return map(FASTA, self.p.raw_dir.contents)

    @property
    def sequences(self):
        """All the sequences from all the files."""
        for fasta in self.raw_files:
            for seq in fasta: yield seq
