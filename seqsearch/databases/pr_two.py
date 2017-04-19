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
class PrTwo(Database):
    """This is the PT2 database.

    https://figshare.com/articles/PR2_rRNA_gene_database/3803709

    To install:
        from seqsearch.databases.pr_two import pr_two
        pr_two.download()
        pr_two.unzip()
        print pr_two.tax_depth_freq

    It will put it in ~/databases/pr_two_11/
    """

    base_url   = "https://ndownloader.figshare.com/articles/3803709/versions/"
    short_name = "pr_two"
    long_name  = 'Protist Ribosomal Reference database (PR2) - SSU rRNA gene database'

    all_paths = """
    /archive.zip
    /pr2_gb203_version_4.5.zip
    /pr2_gb203_version_4.5.fasta
    /pr2_gb203_version_4.5.taxo
    """

    @property
    def rank_names(self):
        """The names of the ranks."""
        return ['Domain',   # 0
                'Kingdom',  # 1
                'Phylum',   # 2
                'Class',    # 3
                'Order',    # 4
                'Family',   # 5
                'Tribe',    # 6
                'Genus',    # 7
                'Species']  # 8

    def __init__(self, version, base_dir=None):
        # Attributes #
        self.version    = version
        self.short_name = self.short_name + "_" + self.version
        # Base directory #
        if base_dir is None: base_dir = home
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # URL #
        self.url = self.base_url + self.version
        # The archive #
        self.dest = self.p.archive
        # The results #
        self.alignment = FilePath(self.base_dir + "pr_two.gb203_v%s.align" % self.version)
        self.taxonomy  = FilePath(self.base_dir + "pr_two.gb203_v%s.tax"   % self.version)
        # The part that mothur will use for naming files #
        self.nickname = "gb203_v%s" % self.version

    def download(self):
        self.dest.directory.create(safe=True)
        self.dest.remove(safe=True)
        print "\nDownloading", self.url
        wget.download(self.url, out=self.dest.path)

    def unzip(self):
        self.dest.unzip_to(self.base_dir, single=False)
        self.p.zip.unzip_to(self.base_dir, single=False)
        self.p.fasta.move_to(self.alignment)
        self.p.taxo.move_to(self.taxonomy)

###############################################################################
pr_two = PrTwo("11")