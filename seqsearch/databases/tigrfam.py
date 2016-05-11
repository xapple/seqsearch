# Built-in modules #

# First party modules #
from seqsearch.databases import Database
from plumbing.cache import property_cached
from fasta import FASTA

###############################################################################
class Tigrfam(Database):
    """The TIGRFAMs is a resource consisting of curated multiple sequence
    alignments, Hidden Markov Models (HMMs) for protein sequence
    classification, and associated information designed to support automated
    annotation of (mostly prokaryotic) proteins.

    http://www.jcvi.org/cgi-bin/tigrfams/index.cgi

    To install:
        from seqsearch.databases.tigrfam import tigrfam
        tigrfam.download()
        tigrfam.unzip()

    It will put it in ~/databases/tigrfam
    """

    all_paths = """
    /raw/
    /unzipped/TIGRFAMs_15.0_HMM.LIB
    /unzipped/TIGRFAMs_15.0_SEED.tar
    /specific/
    """

    short_name = "tigrfam"
    ftp_url    = "ftp.jcvi.org"
    ftp_dir    = "/pub/data/TIGRFAMs/"
    files      = ("TIGRFAMs_15.0_HMM.LIB.gz", "TIGRFAMs_15.0_SEED.tar.gz")

    @property_cached
    def hmm_db(self):
        hmm_db = self.p.hmm
        hmm_db.seq_type = 'hmm_prot'
        return hmm_db

    @property_cached
    def seeds(self):
        seeds = FASTA(self.p.seed)
        return seeds

###############################################################################
tigrfam = Tigrfam("hmm_prot")