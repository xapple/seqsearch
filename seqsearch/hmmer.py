# Futures #
from __future__ import division

# Built-in modules #
import multiprocessing

# Internal modules #
from fasta import FASTA
from plumbing.autopaths import FilePath

# Third party modules #
import sh

###############################################################################
class HMMERquery(object):
    """A Hmmer job."""

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.query)

    def __init__(self, query_path, db_path,
                 seq_type     = 'prot' or 'nucl',     # The seq type of the query_path file
                 params       = None,                 # Add extra params for the command line
                 out_path     = None,                 # Where the results will be dropped
                 executable   = None,                 # If you want a specific binary give the path
                 cpus         = None):                # The number of threads to use
        # Save attributes #
        self.query        = FASTA(query_path)
        self.db           = FilePath(db_path)
        self.seq_type     = seq_type
        self.params       = params if params else {}
        # Output #
        if out_path is None:         self.out_path = self.query.prefix_path + '.hmmout'
        elif out_path.endswith('/'): self.out_path = out_path + self.query.prefix + '.hmmout'
        else:                        self.out_path = out_path
        # Make it a file path #
        self.out_path = FilePath(self.out_path)
        # Executable #
        self.executable = FilePath(executable)
        # Cores to use #
        if cpus is None: self.cpus = multiprocessing.cpu_count()
        else:            self.cpus = cpus

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        else:               cmd = ["hmmersearch"]
        # Essentials #
        cmd += ('-o', '/dev/null',
                '--tblout',  self.p.seq_hits, # parseable table of per-sequence hits
                '--notextw', # unlimited ASCII text output line width
                '--acc',     # prefer accessions over names in output
                '--seed', 1, # set RNG seed to <n>
                self.database,
                self.proteins)
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return map(str, cmd)

    #-------------------------------- RUNNING --------------------------------#
    def run(self):
        """Simply run the HMMER search locally"""
        sh.Command(self.command[0])(self.command[1:])

    @property
    def results(self):
        """Parse the results."""
        for line in self.out_path: yield line.split()