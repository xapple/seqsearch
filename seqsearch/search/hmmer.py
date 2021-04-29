#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import warnings, multiprocessing

# Internal modules #
from seqsearch.databases.pfam    import pfam
from seqsearch.databases.tigrfam import tigrfam

# First party modules #
from fasta import FASTA
from autopaths.file_path import FilePath

# Third party modules #
import sh

# Warnings #
warnings.filterwarnings("ignore", "Bio.SearchIO")
warnings.filterwarnings("ignore", "BiopythonWarning")
from Bio import SearchIO

###############################################################################
class HmmQuery(object):
    """An `hmmsearch` job."""

    short_name = 'hmmsearch'
    long_name  = 'HMMER 3.1b2 (February 2015)'
    executable = 'hmmsearch'
    url        = 'http://hmmer.org/'
    license    = 'GPLv3'
    dependencies = []

    def __nonzero__(self): return bool(self.out_path)

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.query)

    def __init__(self, query_path,                    # The input sequences
                 db_path      = pfam.hmm_db,          # The database to search
                 seq_type     = 'prot' or 'nucl',     # The seq type of the query_path file
                 e_value      = 0.001,                # The search threshold
                 params       = None,                 # Add extra params for the command line
                 out_path     = None,                 # Where the results will be dropped
                 executable   = None,                 # If you want a specific binary give the path
                 cpus         = None):                # The number of threads to use
        # Save attributes #
        self.query      = FASTA(query_path)
        self.db         = FilePath(db_path)
        self.params     = params if params else {}
        self.e_value    = e_value
        self.seq_type   = seq_type
        self.executable = FilePath(executable)
        # Cores to use #
        if cpus is None: self.cpus = min(multiprocessing.cpu_count(), 32)
        else:            self.cpus = cpus
        # Auto detect database short name #
        if db_path == 'pfam':    self.db = pfam.hmm_db
        if db_path == 'tigrfam': self.db = tigrfam.hmm_db
        # Output #
        if out_path is None:
            self.out_path = FilePath(self.query.prefix_path + '.hmmout')
        elif out_path.endswith('/'):
            self.out_path = FilePath(out_path + self.query.prefix + '.hmmout')
        else:
            self.out_path = FilePath(out_path)

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        else:               cmd = ["hmmsearch"]
        # Essentials #
        cmd += ('-o',        '/dev/null',   # direct output to file <f>, not stdout
                '--tblout',  self.out_path, # parsable table of per-sequence hits
                '--seed',    1,             # set RNG seed to <n>
                '--notextw',                # unlimited ASCII text output line width
                '--acc',                    # prefer accessions over names in output
                self.db,
                self.query)
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return map(str, cmd)

    def run(self, cpus=None):
        """Simply run the HMM search locally."""
        # Number of threads #
        if cpus is None: cpus = self.cpus
        # Checks #
        assert self.query.exists
        assert self.db.exists
        # Check if query is not empty #
        if self.query.count_bytes == 0:
            message = "Hmm search on a file with no sequences. File at '%s'"
            warnings.warn(message % self.query, RuntimeWarning)
            return False
        # Do it #
        sh.Command(self.command[0])(['--cpu', str(cpus)] + self.command[1:])

    @property
    def hits(self):
        if not self.out_path:
            raise Exception("You can't access results from HMMER before running the algorithm.")
        return SearchIO.read(self.out_path, 'hmmer3-tab')
