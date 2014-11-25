# Built-in modules #
import math, multiprocessing

# Internal modules #
from seqsearch import SeqSearch
from seqsearch.blast import BLASTquery
from seqsearch.vsearch import VSEARCHquery
from plumbing.cache import property_cached
from fasta.splitable import SplitableFASTA

# Third party modules #
import humanfriendly
from shell_command import shell_output

################################################################################
class ParallelSeqSearch(SeqSearch):
    """The same thing as a SeqSearch but operates by chopping in the input up into
    smaller pieces and running the algorithm on each piece separately, finally joining the outputs.

    You can specify the number of parts, the size in MB or GB that each part should approximately have,
    or even how many sequences should be in each part. Specify only one of the three options.

    In addition, the pieces can be run separately on the local machine, or distributed to different
    compute nodes using the SLURM system.
    """

    def __init__(self, input_fasta, database,
                 num_parts     = None,
                 part_size     = None,
                 seqs_per_part = None,
                 slurm_params  = None,
                 **kwargs):
        # Determine number of parts #
        self.num_parts = None
        # Three possible options #
        if num_parts:
            self.num_parts = num_parts
        if part_size:
            self.bytes_target = humanfriendly.parse_size(part_size)
            self.num_parts = int(math.ceil(input_fasta.count_bytes / self.bytes_target))
        if seqs_per_part:
            self.num_parts = int(math.ceil(input_fasta.count / seqs_per_part))
        # Default case #
        if self.num_parts is None:
            self.num_parts = kwargs.get('num_threads', multiprocessing.cpu_count())
        # In case the user has some special slurm params #
        self.slurm_params = slurm_params
        # Super #
        SeqSearch.__init__(input_fasta, database, **kwargs)

    @property_cached
    def splitable(self):
        """The input fasta file as it is, but with the ability to split it."""
        return SplitableFASTA(self.input_fasta, self.num_parts)

    @property
    def queries(self):
        """A list of all the queries to run."""
        if self.algorithm == 'blast':   return self.blast_queries
        if self.algorithm == 'vsearch': return self.vsearch_queries
        raise NotImplemented(self.algorithm)

    def run(self):
        """Run the search"""
        self.splitable.split()
        for query in self.queries: query.non_block_run()
        for query in self.queries: query.wait()
        self.join_outputs()

    def join_outputs(self):
        """Join the outputs"""
        shell_output('cat %s > %s' % (' '.join(q.out_path for q in self.queries), self.out_path))

    #-------------------------- BLAST IMPLEMTATION -------------------------#
    @property_cached
    def blast_queries(self):
        """Make all BLAST search objects."""
        return [BLASTquery(query_path = p,
                           db_path    = self.database,
                           seq_type   = self.seq_type,
                           params     = self.blast_params,
                           algorithm  = self.select_blast_algo(),
                           version    = "plus",
                           cpus       = 1,
                           num        = p.num) for p in self.splitable.parts]

    #-------------------------- VSEARCH IMPLEMTATION -------------------------#
    @property_cached
    def vsearch_queries(self):
        """Make all VSEARCH search objects."""
        return [VSEARCHquery(p, self.database, self.vsearch_params) for p in self.splitable.parts]