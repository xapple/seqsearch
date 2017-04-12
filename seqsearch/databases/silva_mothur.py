# Built-in modules #
import os, tarfile

# First party modules #
from seqsearch.databases import Database
from fasta import FASTA
from plumbing.autopaths import AutoPaths, FilePath

# Third party modules #
import wget

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class SilvaMothur(Database):
    """This is the SILVA version from mothur.

    https://www.mothur.org/wiki/Silva_reference_files

    To install:
        from seqsearch.databases.silva_mothur import silva_mothur
        silva_mothur.download()
        silva_mothur.unzip()

    It will put it in ~/databases/silva_mothur_xxx/
    """

    base_url   = "https://www.mothur.org/w/images/b/b4/"
    short_name = "silva_mothur"

    all_paths = """
    /silva.tgz
    """

    def __init__(self, version, base_dir=None):
        # Attributes #
        self.version    = version
        self.short_name = self.short_name + "_" + self.version
        # Base directory #
        if base_dir is None: base_dir = home
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # URL #
        self.url = self.base_url + "Silva.nr_v%s.tgz" % self.version
        # The archive #
        self.dest = self.p.tgz
        # The results #
        self.alignment = FilePath(self.base_dir + "silva.nr_v%s.align" % self.version)
        self.taxonomy  = FilePath(self.base_dir + "silva.nr_v%s.tax"   % self.version)
        # The part that mothur will use for naming files #
        self.nickname = "nr_v%s" % self.version

    def download(self):
        self.dest.directory.create(safe=True)
        print "\nDownloading", self.url
        wget.download(self.url, out=self.dest.path)

    def unzip(self):
        # Extract #
        archive = tarfile.open(self.dest, 'r:gz')
        archive.extractall(self.base_dir)
        # Optional testing #
        if self.version == "128": assert self.alignment.md5 == '0b0593da5ee7d1041eb37385f831954e'
        if self.version == "128": assert self.taxonomy.md5  == 'c4d39248146259637e0d0f3ae8bf7e0f'

###############################################################################
silva_mothur = SilvaMothur("128")