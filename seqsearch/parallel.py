# Built-in modules #

# Internal modules #
from seqsearch import SeqSearch
from seqsearch.blast import BLASTquery
from seqsearch.vsearch import VSEARCHquery
from plumbing.cache import property_cached
from fasta.splitable import SplitableFASTA

# Third party modules #
from shell_command import shell_output

################################################################################
class ParallelSeqSearch(SeqSearch):
    """The same thing as a SeqSearch but operates by chopping in the input up into
    smaller pieces and running the algorithm on each piece separately, finally joining the outputs.

    In addition, the pieces can be run separately on the local machine, or distributed to different
    compute nodes using the SLURM system.
    """

    @property_cached
    def splitable(self):
        """The input fasta file, but with the ability to split it."""
        return SplitableFASTA(self.input_fasta, self.num_threads)

    @property
    def queries(self):
        """A list of all the queries to run."""
        if self.algorithm == 'blast':   return self.blast_queries
        if self.algorithm == 'vsearch': return self.vsearch_queries
        raise NotImplemented(self.algorithm)

    @property_cached
    def blast_queries(self):
        """Make all BLAST search objects."""
        # Sequence type #
        if self.seq_type == 'nucl': blast_algo = 'blastn'
        if self.seq_type == 'prot': blast_algo = 'blastp'
        # Make many queries #
        return [BLASTquery(query_path = p,
                           db_path    = self.database,
                           seq_type   = self.seq_type,
                           params     = self.blast_params,
                           algorithm  = blast_algo,
                           version    = "plus",
                           cpus       = 1,
                           num        = p.num) for p in self.splitable.parts]

    @property_cached
    def vsearch_queries(self):
        """Make all VSEARCH search objects."""
        return [VSEARCHquery(p, self.database, self.vsearch_params) for p in self.splitable.parts]

    def run(self):
        """Run the search"""
        self.splitable.split()
        for query in self.queries: query.non_block_run()
        for query in self.queries: query.wait()
        self.join_outputs()

    def join_outputs(self):
        """Join the outputs"""
        shell_output('cat %s > %s' % (' '.join(q.out_path for q in self.queries), self.out_path))