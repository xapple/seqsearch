#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import math, multiprocessing

# Internal modules #
from seqsearch.search import SeqSearch
from seqsearch.search.blast import BLASTquery
from seqsearch.search.vsearch import VSEARCHquery
from plumbing.cache import property_cached
from fasta.splitable import SplitableFASTA

# Third party modules #
from shell_command import shell_output

################################################################################
class ParallelSeqSearch(SeqSearch):
    """
    The same thing as a SeqSearch but operates by chopping in the input up into
    smaller pieces and running the algorithm on each piece separately, finally joining the outputs.

    You can specify the number of parts, the size in MB or GB that each part should approximately have,
    or even how many sequences should be in each part. Specify only one of the three options.

    You can place the pieces in a specific directory.

    In addition, the pieces can be run separately on the local machine, or distributed to different
    compute nodes using the SLURM system.
    """

    def __init__(self, input_fasta, database,
                 num_parts     = None,   # How many fasta pieces should we make
                 part_size     = None,   # What size in MB should a fasta piece be
                 seqs_per_part = None,   # How many sequences in one fasta piece
                 slurm_params  = None,   # Additional parameters for possible SLURM jobs
                 parts_dir     = None,   # If you want a special direcotry for the fasta pieces
                 **kwargs):
        # Determine number of parts #
        self.num_parts = None
        # Three possible options #
        if num_parts:
            self.num_parts = num_parts
        if part_size:
            import humanfriendly
            self.bytes_target = humanfriendly.parse_size(part_size)
            self.num_parts = int(math.ceil(input_fasta.count_bytes / self.bytes_target))
        if seqs_per_part:
            self.num_parts = int(math.ceil(input_fasta.count / seqs_per_part))
        # Default case #
        if self.num_parts is None:
            self.num_parts = kwargs.get('num_threads', min(multiprocessing.cpu_count(), 32))
        # In case the user has some special slurm params #
        self.slurm_params = slurm_params
        # In case the user wants a special parts directory #
        self.parts_dir = parts_dir
        # Super #
        SeqSearch.__init__(self, input_fasta, database, **kwargs)

    @property_cached
    def splitable(self):
        """The input fasta file as it is, but with the ability to split it."""
        return SplitableFASTA(self.input_fasta, self.num_parts, base_dir=self.parts_dir)

    @property
    def queries(self):
        """A list of all the queries to run."""
        if self.algorithm == 'blast':   return self.blast_queries
        if self.algorithm == 'vsearch': return self.vsearch_queries
        raise NotImplemented(self.algorithm)

    def join_outputs(self):
        """Join the outputs"""
        shell_output('cat %s > %s' % (' '.join(q.out_path for q in self.queries), self.out_path))

    #-------------------------------- RUNNING --------------------------------#
    def run_local(self):
        """Run the search locally."""
        # Chop up the FASTA #
        self.splitable.split()
        # Case only one query #
        if len(self.queries) == 1: self.queries[0].run()
        # Case many queries #
        else:
            for query in self.queries: query.non_block_run()
            for query in self.queries: query.wait()
        # Join the results #
        self.join_outputs()

    def run_slurm(self):
        """Run the search via SLURM."""
        self.splitable.split()
        for query in self.queries: query.slurm_job.run()
        for query in self.queries: query.slurm_job.wait()
        self.join_outputs()

    #-------------------------- BLAST IMPLEMENTATION -------------------------#
    @property_cached
    def blast_queries(self):
        """Make all BLAST search objects."""
        return [BLASTquery(query_path   = p,
                           db_path      = self.database,
                           seq_type     = self.seq_type,
                           params       = self.blast_params,
                           algorithm    = self.select_blast_algo(),
                           version      = "plus",
                           cpus         = 1,
                           slurm_params = self.slurm_params,
                           num          = p.num) for p in self.splitable.parts]

    #-------------------------- VSEARCH IMPLEMENTATION -------------------------#
    @property_cached
    def vsearch_queries(self):
        """Make all VSEARCH search objects."""
        return [VSEARCHquery(p, self.database, self.vsearch_params) for p in self.splitable.parts]