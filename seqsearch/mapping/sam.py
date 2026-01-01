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

# Third party modules #
import pysam

###############################################################################
class SamFile(BamOrSamFile):
    """
    Helps you parse a SAM file easily, and offers convenience functions.
    """

    read_mode = 'r'
    format = "SAM"