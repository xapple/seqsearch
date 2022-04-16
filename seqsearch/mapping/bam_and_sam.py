#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import functools, multiprocessing, json

# Internal modules #

# First party modules #
from fasta import PairedFASTQ
from autopaths.file_path import FilePath

# Third party modules #
import pysam

###############################################################################
class BamOrSamFile(FilePath):
    """
    The base class from which both SamFile and BamFile will inherit.
    Conatins method that are relevent for both files
    """

    @functools.cached_property
    def handle(self):
        return pysam.AlignmentFile(self.path, self.read_mode)

    @functools.cached_property
    def __len__(self):
        return self.handle.count()

    def __iter__(self):
        return iter(self.handle)

    @functools.cached_property
    def flagstat(self):
        result = json.loads(pysam.flagstat('-O', 'json', self.path))
        return result['QC-passed reads']

    def didnt_map_to_fastq(self, out_fastq, cpus=None, verbose=True):
        """Extract all the reads that didn't map and save them as a FASTQ."""
        # Message #
        if verbose: print("Extracting non-mapping reads from '%s'." % self)
        # Call method #
        self.extract_to_fastq(out_fastq, 12, None, cpus=cpus)
        # If everything went well, the extra files should be empty #
        assert not out_fastq.singletons
        assert not out_fastq.other_reads
        # So we can remove them #
        out_fastq.singletons.remove()
        out_fastq.other_reads.remove()

    def did_map_to_fastq(self, out_fastq, cpus=None, verbose=True):
        """Extract all the reads that did map and save them as a FASTQ."""
        # Message #
        if verbose: print("Extracting mapping reads from '%s'." % self)
        # Call method #
        self.extract_to_fastq(out_fastq, None, 12, cpus=cpus)
        # If everything went well, the others files should be empty #
        assert not out_fastq.other_reads
        out_fastq.other_reads.remove()

    def extract_to_fastq(self, out_fastq,
                               include = None,
                               exclude = None,
                               cpus    = None,
                               verbose = False):
        """
        Extract a given set of reads and save them as a FASTQ file somewhere
        on the file system.

        Nice introductions and explanations at:
        * https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd
        * https://seqanswers.com/forums/showthread.php?t=17314
        * https://broadinstitute.github.io/picard/explain-flags.html
        """
        # We expect paired FASTQs #
        assert isinstance(out_fastq, PairedFASTQ)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Options #
        options = ['-1', str(out_fastq.fwd),
                   '-2', str(out_fastq.rev),
                   '-s', str(out_fastq.singletons),
                   '-0', str(out_fastq.other_reads),
                   '--threads', str(cpus)]
        # Include #
        if include is not None: options += ['-f', str(include)]
        # Exclude #
        if exclude is not None: options += ['-G', str(exclude)]
        # The input path #
        options += [self.path]
        # Message #
        if verbose:
            msg = ' '.join(options)
            print("Command line is:\n samtools fastq %s \n" % msg)
        # Run our command #
        pysam.fastq(*options)
