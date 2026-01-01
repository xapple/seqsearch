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
from seqsearch.mapping.bam_and_sam import BamOrSamFile

# Third party modules #
import pysam

###############################################################################
class BamFile(BamOrSamFile):
    """
    Helps you parse a BAM file easily, and offers convenience functions.
    """

    read_mode = 'rb'
    format = "BAM"

    def index(self, cpus=None, index=True, verbose=True):
        """Create the corresponding '.bai' file for this BAM file."""
        # Message #
        if verbose: print("Indexing the BAM file '%s'." % self)
        # Number of cores #
        if cpus is None: cpus = min(multiprocessing.cpu_count(), 32)
        # Run the tool #
        return pysam.index('-@', str(cpus), self.path)