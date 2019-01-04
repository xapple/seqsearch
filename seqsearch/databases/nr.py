# Built-in modules #
import urllib
from collections import OrderedDict

# Internal modules #
from seqsearch.blast import BLASTdb
from seqsearch.databases import base_directory

# First party modules #
from fasta import FASTA
from plumbing.autopaths import AutoPaths, FilePath
from plumbing.csv_tables import TSVTable

###############################################################################
class NonRedudant(object):
    """The NR database from NCBI.
    NR contains non-redundant sequences from GenBank translations (i.e. GenPept) together with sequences from other databanks (Refseq, PDB, SwissProt, PIR and PRF)."""

