# Built-in modules #

# Internal modules #
from seqsearch.databases import Database

###############################################################################
class RefSeqBacteriaProtNR(Database):
    """the RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast data base."""

    short_name = "refseq_bact_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/bacteria/"
    pattern    = 'bacteria.nonredundant_protein.*.protein.faa.gz'

###############################################################################
class RefSeqArchaeaProtNR(Database):
    """the RefSeq sequences for only bacteria, only protein, and only the
    non-redundant version. We will download the raw sequences by FTP
    and format them as a blast data base."""

    short_name = "refseq_arch_prot_nr"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/refseq/release/archaea/"
    pattern    = 'archaea.nonredundant_protein.*.protein.faa.gz'

###############################################################################
refseq_bact_prot_nr = RefSeqBacteriaProtNR('prot')
refseq_arch_prot_nr = RefSeqArchaeaProtNR('prot')