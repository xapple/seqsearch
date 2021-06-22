#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing

# First party modules #
from seqsearch.search.blast   import BLASTquery, BLASTdb
from seqsearch.search.vsearch import VSEARCHquery
from seqsearch.search.hmmer   import HmmQuery
from plumbing.cache           import property_cached
from autopaths.file_path      import FilePath

################################################################################
class SeqSearch(object):
    """
    A sequence similarity search.
    Is able to use different algorithms such as BLAST, VSEARCH, HMMER, BLAT etc.
    all through the same interface.

    Input: - Series of sequences in a FASTA file.
           - A database to search against.
           - The type of the sequences ('prot' or 'nucl').
           - The type of algorithm to use. Currently BLAST, VSEARCH or HMMER.
           - Number of threads to use.
           - The desired output path.
           - An extra set of parameters to be given to the search command.
           - The filtering options:
              * BLAST supported:   - e-value
                                   - Maximum targets
                                   - Minimum identity (via manual output format)
                                   - Minimum query coverage (via manual output format)
              * VSEARCH supported: - Maximum targets
              * HMMER supported:   - e-value

    Output: - List of identifiers in the database
              (object with significance value and identity attached)

    Other possible algorithms:

        * https://www.ncbi.nlm.nih.gov/pubmed/11932250
        * https://www.animalgenome.org/bioinfo/resources/manuals/wu-blast/
        * https://github.com/bbuchfink/diamond
    """

    def __repr__(self):
        return '<%s object on %s>' % (self.__class__.__name__, self.input_fasta)

    def __bool__(self):
        return bool(self.query)

    def __init__(self,
                 input_fasta,
                 database,
                 seq_type    = 'prot' or 'nucl',      # What sequence type is the input fasta
                 algorithm   = 'blast' or 'vsearch',  # Which implementation do you want
                 num_threads = None,                  # How many processes to use
                 filtering   = None,                  # The result filtering options
                 out_path    = None,                  # Where the out file will be
                 params      = None,                  # Add extra params for the command line
                 _out        = None,                  # Store the stdout at this path
                 _err        = None):                 # Store the stderr at this path
        # Save attributes #
        self.input_fasta = input_fasta
        self.database    = database
        self.seq_type    = seq_type
        self.algorithm   = algorithm
        self.num_threads = num_threads
        self.filtering   = filtering
        self.out_path    = out_path
        self.params      = params
        self._out        = _out
        self._err        = _err
        # Assign default values #
        self.set_defaults()
        # Validate attributes #
        self.validate()

    def set_defaults(self):
        """
        This method will replace empty attributes with defaults when this is
        needed.
        """
        # In case we got a special object, just use the blast_db attribute #
        if self.algorithm == 'blast' and hasattr(self.database, 'blast_db'):
            self.database = self.database.blast_db
        if self.algorithm == 'vsearch' and hasattr(self.database, 'vsearch_db'):
            self.database = self.database.vsearch_db
        # Otherwise in case we got a path, convert it to a BLASTdb #
        if self.algorithm == 'blast' and not isinstance(self.database, BLASTdb):
            self.database = BLASTdb(self.database)
        # The filtering options #
        if self.filtering is None: self.filtering = {}
        # Output path default value #
        if self.out_path is None:
            self.out_path = self.input_fasta.prefix_path + '.' + \
                            self.algorithm + 'out'
        # Output path setting #
        self.out_path = FilePath(self.out_path)
        # Number of cores default value #
        if self.num_threads is None or self.num_threads is True:
            self.num_threads = min(multiprocessing.cpu_count(), 32)
        # Extra params to be given to the search algorithm #
        if self.params is None: self.params = {}

    def validate(self):
        """
        This method will raise an Exception is any of the arguments passed by
        the user are illegal.
        """
        # The FASTA file has to contain something #
        assert self.input_fasta

    @property
    def query(self):
        """The actual search object with all the relevant parameters."""
        # Pick the right attribute #
        if self.algorithm == 'blast':   return self.blast_query
        if self.algorithm == 'vsearch': return self.vsearch_query
        if self.algorithm == 'hmmer':   return self.hmmer_query
        # Otherwise raise an exception #
        msg = "The algorithm '%s' is not supported."
        raise NotImplemented(msg % self.algorithm)

    #------------------------ Methods that are forwarded ----------------------#
    def run(self):
        """Run the search."""
        return self.query.run()

    def filter(self):
        """Filter the results accordingly."""
        return self.query.filter(self.filtering)

    @property
    def results(self):
        """Parse the results in a meaningful manner."""
        return self.query.results

    #-------------------------- BLAST IMPLEMENTATION -------------------------#
    @property_cached
    def blast_params(self):
        """
        A dictionary of options to pass to the blast executable.
        These params should depend on the filtering options.
        """
        # Make a copy #
        params = self.params.copy()
        # Based on the e-value #
        if 'e_value' in self.filtering:
            params['-evalue'] = self.filtering['e_value']
        # Based on the maximum number of hits #
        if 'max_targets' in self.filtering:
            params['-max_target_seqs'] = self.filtering['max_targets']
        # Return #
        return params

    def select_blast_algo(self):
        """Depends on the query type and the database type."""
        # Pick the right executable #
        db_type = self.database.seq_type
        if self.seq_type == 'nucl' and db_type == 'nucl': return 'blastn'
        if self.seq_type == 'prot' and db_type == 'prot': return 'blastp'
        if self.seq_type == 'nucl' and db_type == 'prot': return 'blastx'
        if self.seq_type == 'prot' and db_type == 'nucl': return 'tblastn'

    @property_cached
    def blast_query(self):
        """Make a BLAST search object."""
        # Make the query object #
        return BLASTquery(query_path = self.input_fasta,
                          db_path    = self.database,
                          seq_type   = self.seq_type,
                          params     = self.blast_params,
                          algorithm  = self.select_blast_algo(),
                          cpus       = self.num_threads,
                          out_path   = self.out_path,
                          _out       = self._out,
                          _err       = self._err)

    #------------------------- VSEARCH IMPLEMENTATION ------------------------#
    @property_cached
    def vsearch_params(self):
        """
        A dictionary of options to pass to the vsearch executable.
        These params should depend on the filtering options.
        """
        # Make a copy #
        params = self.params.copy()
        # Based on the maximum number of hits #
        if 'max_targets' in self.filtering:
            params['-maxaccepts'] = self.filtering['max_targets']
        # Return #
        return params

    @property_cached
    def vsearch_query(self):
        """Make a VSEARCH search object."""
        return VSEARCHquery(query_path = self.input_fasta,
                            db_path    = self.database,
                            seq_type   = self.seq_type,
                            params     = self.vsearch_params,
                            algorithm  = "usearch_global",
                            cpus       = self.num_threads,
                            out_path   = self.out_path,
                            _out       = self._out,
                            _err       = self._err)

    #-------------------------- HMMER IMPLEMENTATION -------------------------#
    @property_cached
    def hmmer_params(self):
        params = self.params.copy()
        if 'e_value' in self.filtering: params['-E'] = self.filtering['e_value']
        return params

    @property_cached
    def hmmer_query(self):
        """Make a HMMER search object."""
        return HmmQuery(query_path = self.input_fasta,
                        db_path    = self.database,
                        seq_type   = self.seq_type,
                        params     = self.hmmer_params,
                        cpus       = self.num_threads,
                        out_path   = self.out_path)
