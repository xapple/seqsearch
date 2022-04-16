#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, tarfile, shutil

# First party modules #
from seqsearch.databases import DatabaseHTTP
from fasta import FASTA
from plumbing.timer import Timer

# Third party modules #
import pandas
from tqdm import tqdm

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class HumGut(DatabaseHTTP):
    """
    Quote from the paper:

    "We screened > 5,700 healthy human gut metagenomes for the containment of
    > 490,000 publicly available prokaryotic genomes sourced from RefSeq and
    the recently announced UHGG collection. This resulted in a pool of
    > 381,000 genomes that were subsequently scored and ranked based on their
    prevalence in the healthy human metagenomes. The genomes were then
    clustered at a 97.5% sequence identity resolution, and cluster
    representatives (30,691 in total) were retained to comprise the HumGut
    collection."

    The publication is here:

    * https://doi.org/10.1186/s40168-021-01114-w

    The download website is here:

    * http://arken.nmbu.no/~larssn/humgut/

    To install:

        >>> from seqsearch.databases.humgut import humgut
        >>> humgut.download()
        >>> humgut.untargz()
        >>> humgut.get_95_cluster()
        >>> humgut.autopaths.tar.remove()
        >>> humgut.make_bwa_database()

    It will place the resulting files in "~/databases/human/".
    """

    tag        = "humgut"
    short_name = "humgut"
    long_name  = 'HumGut: a comprehensive human gut prokaryotic genomes' \
                 ' collection'

    base_url = "http://arken.nmbu.no/~larssn/humgut/"

    files = ['HumGut.tar.gz',
             'HumGut.tsv',
             'ncbi_names.dmp',
             'ncbi_nodes.dmp']

    all_paths = """
    /raw/HumGut.tar.gz
    /raw/HumGut.tsv
    /raw/ncbi_names.dmp
    /raw/ncbi_nodes.dmp
    /cluster95/humgut95.fasta.gz
    /bwa_db/
    """

    def get_95_cluster(self, verbose=True):
        """
        See the documentation on GitHub here:
        * https://github.com/larssnip/HumGut#the-humgut-library
        """
        # Parse the TSV metadata with pandas #
        df = pandas.read_csv(self.autopaths.HumGut_tsv, sep='\t')
        # Keep the first representative of every cluster #
        df = df.drop_duplicates(subset = 'cluster95')
        # Get the list of genome names to retrieve #
        genomes_to_get = list(df['genome_file'])
        # Function #
        def find_entry(name, all_entries):
            for entry in all_entries:
                if entry.name == name:
                    return entry
            raise Exception("Entry '%s' not found." % name)
        # Open the tar file with all genomes for reading.
        # Fetch every one of the separate FASTA files (one per genome) from
        # within the archive and concatenate them all into one big FASTA file.
        with tarfile.open(self.autopaths.tar, "r:gz") as tar:
            # Message #
            msg = "Reading file list from '%s'..."
            print(msg % self.autopaths.tar)
            # Get all the member files #
            members = tar.getmembers()
            # Message #
            msg = "Extracting %i genomes..."
            print(msg % len(genomes_to_get))
            # Iterate #
            with open(self.autopaths.fasta, "wb") as out_file:
                for genome in tqdm(genomes_to_get):
                    info = find_entry('fna/' + genome, members)
                    handle = tar.extractfile(info)
                    shutil.copyfileobj(handle, out_file)
                    break

    @property
    def bwa_index(self):
        return self.autopaths.bwa_db_dir + 'humgut95'

    def make_bwa_database(self, verbose=True, print_time=True):
        """
        Using the 95% clustered FASTA file, create a BWA compatible database.
        On a typical single threaded Intel process this takes about:
            [main] Real time: 22335.125 sec; CPU: 22106.508 sec
        """
        # The big fasta with all the genomes #
        fasta = FASTA(self.autopaths.fasta)
        # Create a timer #
        if print_time:
            timer = Timer()
            timer.print_start()
        # Make a BWA index with the 'bwtsw' algorithm #
        fasta.index_bwa(self.bwa_index, verbose=verbose)
        # End message #
        if print_time:
            timer.print_end()
            timer.print_total_elapsed()

###############################################################################
# Create a singleton #
humgut = HumGut()