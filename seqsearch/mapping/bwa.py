#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import functools, multiprocessing, sys

# Internal modules #
from seqsearch.mapping.sam import SamFile, BamFile

# First party modules #
from plumbing.check_cmd_found import check_cmd
from fasta import PairedFASTQ

# Third party modules #
import sh

###############################################################################
class MapBWA:
    """
    Takes care of running the `bwa` program to map reads from a FASTQ file
    pair to a reference (e.g. a genome).

    Installing:
    * If you are on Ubuntu you can just type "apt install bwa".
    * If you are on macOS you can just type "brew install bwa".
    """

    def __repr__(self):
        msg = '<%s object on "%s">'
        return msg % (self.__class__.__name__, self.source)

    def __init__(self, source, ref_db, sam_out):
        # Source is a FASTQ file that will contain the reads to map #
        assert isinstance(source, PairedFASTQ)
        self.source = source
        # The ref is a BWA formatted database that will contain the reference #
        self.ref_db = ref_db
        # The destination is a SAM file that will contain the alignments #
        self.sam_out = SamFile(sam_out)
        # For convenience, we will provide the path of the matching BAM file
        # However, this one needs to be created separately at a later time.
        self.bam = BamFile(sam_out.replace_extension('bam'))

    #------------------------------ Running ----------------------------------#
    def __call__(self, cpus=None, verbose=True, stderr=False, timer=False):
        """Run BWA MEM, and produce a SAM file."""
        # Message #
        if verbose:
            msg = "\n\nRunning `bwa mem` from '%s' onto '%s'.\n"
            print(msg % (self.source, self.ref_db.bwa_index))
        # Timer #
        if timer:
            from plumbing.timer import Timer
            timer = Timer()
        # Check it is installed #
        check_cmd('bwa')
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Get the command #
        bwa = sh.Command("bwa")
        # You can't specify an output file, it always uses the stdout #
        options = {'t':    cpus,
                   '_out': self.sam_out}
        # Output the status messages #
        if stderr: options['_err'] = sys.stderr
        # Run it #
        bwa('mem',
            self.ref_db.bwa_index,
            self.source.fwd,
            self.source.rev,
            **options)
        # Elapsed time #
        if timer: timer.print_total_elapsed(msg='BWA: ')
        # Return #
        return self.sam_out

    #------------------------------- Results ---------------------------------#
    def __bool__(self):
        """
        Return True if the `bwa` software was run already and the
        results are stored on the filesystem. Return False if it was not yet
        run.
        """
        return self.sam_out.exists

    @functools.cached_property
    def results(self):
        # Check it was run #
        if not self:
            msg = "You can't access results from `bwa` " \
                  "before running the tool."
            raise Exception(msg)
        # Return the results #
        return self.sam_out
