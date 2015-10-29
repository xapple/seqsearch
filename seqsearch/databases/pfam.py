# Built-in modules #
import os

# Internal modules #
from seqsearch.databases import Database

# Constants #
home = os.environ['HOME'] + '/'

###############################################################################
class Pfam(Database):
    """The Pfam database is a large collection of protein families,
    each represented by multiple sequence alignments and HMMs.
    http://pfam.xfam.org
    Pfam-A is the manually curated portion of the database that
    contains over 16,000 entries. An automatically generated supplement
    is provided called Pfam-B. Pfam-B was discontinued."""

    short_name = "pfam"
    ftp_url    = "ftp.ebi.ac.uk"
    ftp_dir    = "/pub/databases/Pfam/releases/Pfam28.0/"
    pattern    = 'Pfam-A.hmm.gz'


    @property
    def hmm_db(self):
        pass

    def unzip(self):
        """Unzip them"""
        for f in self.raw_files: f.ungzip_to(self.p.unzipped_dir + f.prefix)

###############################################################################
pfam = Pfam("hmm")