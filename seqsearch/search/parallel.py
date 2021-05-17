#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import math, subprocess, multiprocessing

# First party modules #
from seqsearch.search         import SeqSearch
from seqsearch.search.blast   import BLASTquery
from seqsearch.search.vsearch import VSEARCHquery
from plumbing.cache           import property_cached
from fasta.splitable          import SplitableFASTA

################################################################################
class ParallelSeqSearch(SeqSearch):
    """
    The same thing as a SeqSearch but operates by chopping in the input up into
    smaller pieces and running the algorithm on each piece separately, finally
    joining the outputs.

    You can specify the number of parts, the size in MB or GB that each part
    should approximately have, or even how many sequences should be in each
    part. Specify only one of the three options.

    You can place the pieces in a specific directory.
    """

    def __init__(self,
                 input_fasta,
                 database,
                 num_parts     = None,   # How many fasta pieces should we make
                 part_size     = None,   # What size in MB should a fasta piece be
                 seqs_per_part = None,   # How many sequences in one fasta piece
                 parts_dir     = None,   # If you want a special directory for the fasta pieces
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
            default = min(multiprocessing.cpu_count(), 32)
            self.num_parts = kwargs.get('num_threads', default)
            if self.num_parts is True: self.num_parts = default
        # In case the user wants a special parts directory #
        self.parts_dir = parts_dir
        # Super #
        SeqSearch.__init__(self, input_fasta, database, **kwargs)

    @property_cached
    def splitable(self):
        """The input fasta file as it is, but with the ability to split it."""
        return SplitableFASTA(self.input_fasta,
                              self.num_parts,
                              base_dir = self.parts_dir)

    @property
    def queries(self):
        """A list of all the queries to run."""
        if self.algorithm == 'blast':   return self.blast_queries
        if self.algorithm == 'vsearch': return self.vsearch_queries
        raise NotImplemented(self.algorithm)

    def join_outputs(self):
        """Join the outputs."""
        all_files = ' '.join(q.out_path for q in self.queries)
        subprocess.check_call('cat %s > %s' % (all_files, self.out_path))

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

    #-------------------------- BLAST IMPLEMENTATION -------------------------#
    @property_cached
    def blast_queries(self):
        """Make all BLAST search objects."""
        # Select the right BLAST #
        blast_algo = self.select_blast_algo()
        # Create a list #
        return [BLASTquery(query_path   = p,
                           db_path      = self.database,
                           seq_type     = self.seq_type,
                           params       = self.blast_params,
                           algorithm    = blast_algo,
                           cpus         = 1,
                           num          = p.num) for p in self.splitable.parts]

    #-------------------------- VSEARCH IMPLEMENTATION -------------------------#
    @property_cached
    def vsearch_queries(self):
        """Make all VSEARCH search objects."""
        return [VSEARCHquery(p, self.database, self.vsearch_params)
                for p in self.splitable.parts]