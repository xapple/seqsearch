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
    is provided called Pfam-B. Pfam-B was discontinued.

    To install:
        from seqsearch.databases.pfam import pfam
        pfam.download()
        pfam.unzip()
    """

    short_name = "pfam"
    ftp_url    = "ftp.ebi.ac.uk"
    ftp_dir    = "/pub/databases/Pfam/releases/Pfam28.0/"
    pattern    = 'Pfam-A.hmm.gz'

    @property
    def hmm_db(self):
        return self.p.unzipped_dir.contents.next()

###############################################################################
pfam = Pfam("hmm")