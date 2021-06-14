#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing, threading

# First party modules #
from autopaths.file_path import FilePath
from fasta import FASTA

###############################################################################
class CoreSearch(object):
    """
    A class to inherit from.
    Contains methods that are common to all search algorithms implementation.
    Currently: BLASTquery and VSEARCHquery inherit from this.
    """

    extension = 'out'

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.query)

    def __bool__(self):
        return bool(self.out_path)

    def __init__(self,
                 query_path,
                 db_path,
                 seq_type     = 'prot' or 'nucl',     # The seq type of the query_path file
                 params       = None,                 # Add extra params for the command line
                 algorithm    = "blastn" or "blastp", # Will be auto-determined with seq_type
                 out_path     = None,                 # Where the results will be dropped
                 executable   = None,                 # If you want a specific binary give the path
                 cpus         = None,                 # The number of threads to use
                 num          = None,                 # When parallelized, the number of this thread
                 _out         = None,                 # Store the stdout at this path
                 _err         = None):                # Store the stderr at this path
        # Main input #
        self.query = FASTA(query_path)
        # The database to search against #
        self.db = FilePath(db_path)
        # Other attributes #
        self.seq_type     = seq_type
        self.algorithm    = algorithm
        self.num          = num
        self.params       = params if params else {}
        # The standard output and error #
        self._out         = _out
        self._err         = _err
        # Output defaults #
        if out_path is None:
            self.out_path = self.query.prefix_path + self.extension
        elif out_path.endswith('/'):
            self.out_path = out_path + self.query.prefix + self.extension
        else:
            self.out_path = out_path
        # Make it a file path #
        self.out_path = FilePath(self.out_path)
        # Executable #
        self.executable = FilePath(executable)
        # Cores to use #
        if cpus is None: self.cpus = min(multiprocessing.cpu_count(), 32)
        else:            self.cpus = cpus
        # Save the output somewhere #
        if self._out is True:
            self._out = self.out_path + '.stdout'
        if self._err is True:
            self._err = self.out_path + '.stderr'

    #-------------------------------- RUNNING --------------------------------#
    def non_block_run(self):
        """Special method to run the query in a thread without blocking."""
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True # So that they die when we die
        self.thread.start()

    def wait(self):
        """
        If you have run the query in a non-blocking way, call this method to pause
        until the query is finished.
        """
        try:
            # We set a large timeout so that we can Ctrl-C the process
            self.thread.join(999999999)
        except KeyboardInterrupt:
            print("Stopped waiting on BLAST thread number %i" % self.num)

