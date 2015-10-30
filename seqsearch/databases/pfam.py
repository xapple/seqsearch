# Built-in modules #
import os, sh

# First party modules #
from seqsearch.databases import Database
from plumbing.autopaths import DirectoryPath, AutoPaths
from plumbing.cache import property_cached
from fasta import FASTA

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
    /unzipped/Pfam-A.hmm
    /unzipped/Pfam-A.fasta
    /specific/
    """

    short_name = "pfam"
    ftp_url    = "ftp.ebi.ac.uk"
    ftp_dir    = "/pub/databases/Pfam/current_release/"
    files      = ("Pfam-A.hmm.gz", "Pfam-A.fasta.gz")

    @property_cached
    def hmm_db(self):
        hmm_db = self.p.hmm
        hmm_db.seqtype = 'hmm_prot'
        return hmm_db

    @property_cached
    def fasta(self):
        fasta = FASTA(self.p.fasta)
        return fasta

###############################################################################
pfam = Pfam("hmm")

###############################################################################
class SpecificFamily(object):
    """When you are interested in having an HMM 'database' with only
    one specific Pfam in it."""

    all_paths = """
    /model.hmm
    /proteins.fasta
    """

    def __init__(self, fam_name):
        self.fam_name = fam_name
        self.base_dir = DirectoryPath(pfam.p.specific_dir + self.fam_name)
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def hmm_db(self):
        hmm_db = self.p.model
        hmm_db.seqtype = 'hmm_prot'
        if not hmm_db.exists:
            print sh.hmmfetch('-o', hmm_db, pfam.hmm_db, self.fam_name)
            assert hmm_db
        return hmm_db

    @property_cached
    def fasta(self):
        fasta = FASTA(self.p.proteins)
        if not fasta.exists:
            fasta.create()
            for seq in pfam.fasta:
                if self.fam_name in seq.description: fasta.add_seq(seq)
            fasta.close()
            assert fasta
        return fasta
