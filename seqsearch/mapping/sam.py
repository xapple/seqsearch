#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import multiprocessing

# Internal modules #
from seqsearch.mapping.bam import BamFile
from seqsearch.mapping.bam_and_sam import BamOrSamFile

# First party modules #

# Third party modules #
import pysam

###############################################################################
class SamFile(BamOrSamFile):
    """
    Helps you parse a SAM file easily, and offers convenience functions.
    """

    read_mode = 'r'

    def collate_to_bam(self, out_bam, cpus=None, verbose=True):
        """
        Collate the current SAM file and produce an uncompressed BAM.
        Subsequently, running `samtols flagstats` produces the same results
        between the original SAM and the collated BAM.
        """
        # Message #
        if verbose: print("Collating the SAM file '%s'." % self)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Make sure it is a BAM object #
        out_bam = BamFile(out_bam)
        # Run #
        pysam.collate('-o',           str(out_bam),
                      '--output-fmt', 'bam',
                      '--threads',    str(cpus),
                      '-u',           self.path)
        # Return #
        return out_bam

    def sort_to_bam(self, out_bam, cpus=None, verbose=True):
        """Sort the current SAM file and produce a compressed BAM file."""
        # Message #
        if verbose: print("Sorting the SAM file '%s'." % self)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Make sure it is a BAM object #
        out_bam = BamFile(out_bam)
        # Run the sorting #
        pysam.sort('-o',           str(out_bam),
                   '--output-fmt', 'bam',
                   '--threads',    str(cpus),
                   self.path)
        # Return #
        return out_bam