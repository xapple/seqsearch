b'This module needs Python 2.7.x'

# Special variables #
__version__ = '0.0.1'

# Built-in modules #
import multiprocessing

# Internal modules #
from seqsearch.blast import BLASTquery
from seqsearch.vsearch import VSEARCHquery
from plumbing.cache import property_cached
from plumbing.autopaths import FilePath

# Third party modules #

################################################################################
class SeqSearch(object):
    """A sequence similarity search. Could use different algorithms such
    as BLAST, USEARCH, BLAT etc.

    Input: - List of sequences in a FASTA file
           - The type of the sequences
           - A database to search against
           - The type of algorithm to use
           - Number of threads to use
           - The desired output path
           - The filtering options:
             * BLAST supported:   - Minimum identity
                                  - E value
                                  - Maximum targets
                                  - Minimum query coverage (via manual output format)
             * USEARCH supported: - ?
    Output: - Sorted list of identifiers in the database (object with significance value and identity attached)
    """

    def __init__(self, input_fasta, seq_type, database,
                 algorithm   = 'blast',
                 num_threads = None,
                 filtering   = None,
                 out_path    = None,
                 params      = None):
        # Base parameters #
        self.input_fasta = input_fasta
        self.seq_type = seq_type
        self.database = database
        # Optional #
        self.algorithm = algorithm
        # The filtering options #
        if filtering is None: self.filtering = {}
        else: self.filtering = filtering
        # Number of cores to use #
        if num_threads is None: self.num_threads = multiprocessing.cpu_count()
        else: self.num_threads = num_threads
        # Output path #
        if out_path is None: self.out_path = FilePath(self.input_fasta.prefix_path + '.' + algorithm + 'out')
        else: self.out_path = FilePath(out_path)
        # Extra params to be past to the search algorithm #
        if num_threads is None: self.params = {}
        else: self.params = params

    @property
    def query(self):
        """The similarity search object with all the relevant parameters."""
        if self.algorithm == 'blast':   return self.blast_query
        if self.algorithm == 'usearch': return self.vsearch_query
        raise NotImplemented(self.algorithm)

    @property_cached
    def blast_params(self):
        """A dictionary of options to pass to the blast executable.
        The params should depend on the filtering options."""
        # Initialize #
        params = self.params.copy()
        # Defaults #
        params['-num_threads'] = self.num_threads
        # Conditionals #
        if 'e_value'      in self.filtering: params['-evalue']          = self.filtering['e_value']
        if 'max_targets'  in self.filtering: params['-max_target_seqs'] = self.filtering['max_targets']
        # Depends on sequence type #
        if self.seq_type == 'nucl':
            if 'min_identity' in self.filtering: params['-perc_identity'] = self.filtering['min_identity'] * 100
        # Return #
        return params

    @property_cached
    def blast_query(self):
        """Make a BLAST search object."""
        # Sequence type #
        if self.seq_type == 'nucl': blast_algo = 'blastn'
        if self.seq_type == 'prot': blast_algo = 'blastp'
        # The query object #
        return BLASTquery(query_path = self.input_fasta,
                          db_path    = self.database,
                          params     = self.blast_params,
                          algorithm  = blast_algo,
                          version    = "plus",
                          out_path   = self.out_path)

    @property_cached
    def vsearch_query(self):
        """Make a VSEARCH search object."""
        params = {}
        query = VSEARCHquery(self.input_fasta, self.database, params)
        return query

    def run(self):
        """Run the search"""
        return self.query.run()

    def filter(self):
        """Filter the results accordingly"""
        return self.query.filter(self.filtering)

    @property
    def results(self):
        """Parse the results."""
        return self.query.results