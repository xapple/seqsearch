# Built-in modules #
import os

# First party modules #
from seqsearch.databases import Database
from plumbing.cache import property_cached
from fasta import FASTA
from plumbing.autopaths import AutoPaths, FilePath

# Third party modules #
import wget

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class Silva(Database):
    """SILVA provides comprehensive, quality checked and regularly updated datasets of aligned small (16S/18S, SSU) and large subunit (23S/28S, LSU) ribosomal RNA (rRNA) sequences for all three domains of life (Bacteria, Archaea and Eukarya). SILVA are the official databases of the software package ARB.

    https://www.arb-silva.de

    To install:
        from seqsearch.databases.silva import Silva
        tigrfam.download()
        tigrfam.unzip()

    It will put it in ~/databases/silva_nnn/
    """

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
        # The database #
        self.name = "SILVA_%s.1_SSURef_Nr99_tax_silva.fasta.gz" % self.version
        self.url  = "release_%s_1/Exports/"  % self.version
        self.dest = FilePath(self.base_dir + self.name)

    def download(self):
        wget.download(self.base_url + self.url + self.name, self.dest)

###############################################################################
silva = Silva("123", "nucl")