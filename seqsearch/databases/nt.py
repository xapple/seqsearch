#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# Internal modules #
from seqsearch.databases import Database
from seqsearch.search    import SeqSearch

# First party modules #
from autopaths.tmp_path import new_temp_dir
from fasta import FASTA

###############################################################################
class NucleotideDatabase(Database):
    """
    The Nucleotide database is a collection of sequences from several sources,
     including GenBank, RefSeq, TPA and PDB.

     To install:

        from seqsearch.databases.nt import nt
        nt.download()
        nt.untargz()
        nt.test()

    It will put it in ~/databases/nt
    """

    short_name = "nt"
    long_name  = "The Nucleotide database (NCBI)"
    ftp_url    = "ftp.ncbi.nlm.nih.gov"
    ftp_dir    = "/blast/db/"
    pattern    = 'nt.*.tar.gz'

    def test(self):
        """Search one sequence, and see if it works."""
        # New directory #
        directory = new_temp_dir()
        # A randomly chosen sequence (Homo sapiens mRNA for prepro cortistatin) #
        seq = """ACAAGATGCCATTGTCCCCCGGCCTCCTGCTGCTGCTGCTCTCCGGGGCCACGGCCACCGCTGCCCTGCC
        CCTGGAGGGTGGCCCCACCGGCCGAGACAGCGAGCATATGCAGGAAGCGGCAGGAATAAGGAAAAGCAGC
        CTCCTGACTTTCCTCGCTTGGTGGTTTGAGTGGACCTCCCAGGCCAGTGCCGGGCCCCTCATAGGAGAGG
        AAGCTCGGGAGGTGGCCAGGCGGCAGGAAGGCGCACCCCCCCAGCAATCCGCGCGCCGGGACAGAATGCC
        CTGCAGGAACTTCTTCTGGAAGACCTTCTCCTCCTGCAAATAAAACCTCACCCATGAATGCTCACGCAAG
        TTTAATTACAGACCTGAA"""
        seq = seq.replace('\n','')
        seq = seq.replace(' ','')
        # Make input #
        input_fasta = FASTA(directory + 'input.fasta')
        input_fasta.create()
        input_fasta.add_str(seq, "My test sequence")
        input_fasta.close()
        # Make output #
        out_path = directory + 'output.blast'
        # Make extras parameters #
        params = {'-outfmt': 0,
                  '-evalue': 1e-5,
                  '-perc_identity': 99}
        # Make the search #
        search = SeqSearch(input_fasta,
                           self.blast_db,
                           'nucl',
                           'blast',
                           num_threads = 1,
                           out_path    = out_path,
                           params      = params)
        # Run it #
        search.run()
        # Print result #
        print("Success", directory)

###############################################################################
nt = NucleotideDatabase("nucl")