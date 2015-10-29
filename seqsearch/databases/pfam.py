# Built-in modules #
import os, sh

# Internal modules #
from seqsearch.databases import Database
from plumbing.autopaths import DirectoryPath, FilePath, AutoPaths

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

    all_paths = """
    /raw/
    /unzipped/
    /specific/
    """

    short_name = "pfam"
    ftp_url    = "ftp.ebi.ac.uk"
    ftp_dir    = "/pub/databases/Pfam/current_release/"
    pattern    = 'Pfam-A.hmm.gz'

    @property
    def hmm_db(self):
        return self.p.unzipped_dir.contents.next()

###############################################################################
pfam = Pfam("hmm")

###############################################################################
class SpecificFamily(Database):
    """When you are interested in having an HMM 'database' with only
    one specific Pfam in it.
    """

    all_paths = """
    /model.hmm
    """

    def __init__(self, fam_name):
        self.fam_name = fam_name
        self.base_dir = DirectoryPath(pfam.p.specific_dir + self.fam_name)
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property
    def hmm_db(self):
        hmm_db = self.p.model
        if not hmm_db.exists:
            print sh.hmmfetch('-o', hmm_db, pfam.hmm_db, self.fam_name)
            assert hmm_db
        return hmm_db
