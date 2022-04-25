#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, re, tarfile, shutil, functools

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
class HumGutGenome:
    """A genome found in the database, with associated counts."""

    def __init__(self, tax_id, count):
        """Use the HumGut genome ID to create the instance."""
        # Save attributes #
        self.tax_id = tax_id
        self.count  = count

    def __repr__(self):
        return '<%s object ID %s>' % (self.__class__.__name__, self.tax_id)

    @functools.cached_property
    def metadata(self):
        """
        An example output is:

            {'HumGut_name': 'HumGut_20705',
             'cluster975': 20705,
             'cluster95': 3214,
             'gtdbtk_tax_id': 4030631,
             'gtdbtk_organism_name': 's__Enterococcus_D casseliflavus',
             'gtdbtk_taxonomy': 'd__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus_D;s__Enterococcus_D casseliflavus',
             'ncbi_tax_id': 1218087,
             'ncbi_organism_name': 'Enterococcus casseliflavus NBRC 100478',
             'ncbi_rank': 'strain',
             'prevalence_score': 0.7849960826259196,
             'metagenomes_present': 22,
             'completeness': 99.24528301886792,
             'contamination': 0.389531345100426,
             'GC': 0.4235326316891364,
             'genome_size': 3668336,
             'source': 'RefSeq',
             'genome_type': 'Complete Genome',
             'cluster975_size': 68,
             'cluster95_size': 125,
             'genome_file': 'GCF_003641225.1_ASM364122v1_genomic.fna.gz',
             'ftp_download': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/641/225/GCF_003641225.1_ASM364122v1/GCF_003641225.1_ASM364122v1_genomic.fna.gz'}
        """
        return humgut.id_to_metadata(self.tax_id)

    @functools.cached_property
    def tax(self):
        """
        Parse the gtdbtk_taxonomy string into a list.
        A typical output is the following:

            ['Bacteria',
             'Firmicutes',
             'Bacilli',
             'Lactobacillales',
             'Enterococcaceae',
             'Enterococcus_D']
        """
        # Parse the string 'd__Bacteria;p__Firmicutes;...' #
        pattern = '__(.+?);'
        tax = re.findall(pattern, self.metadata['gtdbtk_taxonomy'])
        # If the classification doesn't go all the way down we will add
        # 'Unclassified' until reaching the species level
        count_ranks = len(humgut.rank_names)
        tax += ['Unclassified'] * (count_ranks - len(tax))
        # Return tax #
        return tax

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

    #------------------------------ Properties -------------------------------#
    @functools.cached_property
    def metadata(self):
        """Parse the TSV metadata with pandas"""
        # Load the file in memory #
        df = pandas.read_csv(self.autopaths.HumGut_tsv, sep='\t')
        # Change index #
        df = df.set_index('HumGut_tax_id')
        # Return #
        return df

    @property
    def bwa_index(self):
        return self.autopaths.bwa_db_dir + 'humgut95'

    #--------------------------- Extra information ----------------------------#
    @property
    def rank_names(self):
        return ['Domain',   # 1 (This is Bacteria, Archaea or Eucarya)
                'Phylum',   # 2 (This is for instance 'Firmicutes')
                'Class',    # 3
                'Order',    # 4
                'Family',   # 5
                'Genus',    # 6
                'Species']  # 7

    #------------------------------- Methods ---------------------------------#
    def id_to_metadata(self, tax_id):
        """
        Like a dictionary for quick look up of a genome based on its
        taxonomy ID.
        """
        return self.metadata.loc[tax_id]

    def get_95_cluster(self):
        """
        Will write a large compressed FASTA file will all genomes of
        interest concatenated together. To do this, we will pick single files
        out of the large TAR archive provided and append them one by one.
        See the documentation on GitHub here:
        * https://github.com/larssnip/HumGut#the-humgut-library
        """
        # Keep the first representative of every cluster #
        df = self.metadata.drop_duplicates(subset = 'cluster95')
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