#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os

# First party modules #
from seqsearch.databases import Database
from fasta import FASTA
from autopaths.auto_paths import AutoPaths
from autopaths.file_path import FilePath

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class Foraminifera(Database):
    """
    This is a custom database containing exclusively Foraminifera sequences.

    https://genev.unige.ch/research/laboratory/Jan-Pawlowski

    You should place the file "foram_db_cor.fasta" in:  ~/databases/foraminifera/
    Then you can run this:

            from seqsearch.databases.foraminifera import foraminifera
            foraminifera.process()
            print foraminifera.tax_depth_freq

    """

    short_name = "foraminifera"
    long_name  = 'The custom made Foraminifera database as received by email on 7th April 2017'

    all_paths = """
    /foram_db_cor.fasta
    /foram_mothur.fasta
    /foram_mothur.tax
    """

    @property
    def rank_names(self):
        """The names of the ranks. Total 9 ranks."""
        return ['Domain',   # 0
                'Kingdom',  # 1
                'Phylum',   # 2
                'Class',    # 3
                'Order',    # 4
                'Family',   # 5
                'Tribe',    # 6
                'Genus',    # 7
                'Species']  # 8

    def __init__(self, base_dir=None):
        # Base directory #
        if base_dir is None: base_dir = home
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # The results #
        self.alignment = FASTA(self.p.mothur_fasta)
        self.taxonomy  = FilePath(self.p.mothur_tax)
        # The part that mothur will use for naming files #
        self.nickname = "foram_mothur"

    def process(self):
        # The file that was received by email without documentation #
        raw = FASTA(self.p.cor)
        # Open files #
        self.alignment.create()
        self.taxonomy.create()
        # Loop #
        for seq in raw:
            # Parse #
            name = seq.id[11:].split('|')
            num  = name.pop(0)
            # Check #
            for x in name: assert ';' not in x
            for x in name: assert '\t' not in x
            # Make ranks #
            ranks = ['Eukaryota'                       , # 0 Domain
                     'Rhizaria'                        , # 1 Kingdom
                     'Foraminifera'                    , # 2 Phylum
                     name[0]                           , # 3 Class
                     name[1]                           , # 4 Order
                     name[2]                           , # 5 Family
                     name[3]                           , # 6 Tribe
                     name[4]                           , # 7 Genus
                     name[5]]                            # 8 Species
            # The taxonomy string #
            tax_line = ';'.join(ranks)
            # Add sequence to the new fasta file #
            self.alignment.add_str(str(seq.seq), name="foram" + num)
            # Add the taxonomy to the tax file #
            self.taxonomy.add_str("foram" + num + '\t' + tax_line + '\n')
        # Close files #
        self.alignment.close()
        self.taxonomy.close()

###############################################################################
foraminifera = Foraminifera()