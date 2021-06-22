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

# Internal modules #
from Bio import SearchIO
from autopaths.file_path import FilePath
from seqsearch.search.core import CoreSearch

# Module for launching shell commands #
if os.environ.get('INSIDE_PYCHARM'): import pbs3 as sh
else:                                import sh

###############################################################################
class VSEARCHquery(CoreSearch):
    """
    A vsearch job.

    Example commands:

        vsearch --usearch_global queries.fas --db queries.fas --id 0.6
                --blast6out results.blast6 --maxaccepts 0 --maxrejects 0

        ./vsearch --usearch_global queries.fsa --db database.fsa --id 0.9
                  --alnout alnout.txt

    Some interesting options form the manual follow:

    --maxaccepts <positive integer>
    > Maximum number of hits to accept before stopping the search. The default
    value is 1. This option works in pair with --maxrejects. The search process
    sorts target sequences by decreasing number of k-mers they have in common
    with the query sequence, using that information as a proxy for sequence
    similarity. After pairwise alignments, if the first target sequence passes
    the acceptation criteria, it is accepted as best hit and the search process
    stops for that query. If --maxaccepts is set to a higher value, more hits
    are accepted.
    """

    extension = 'vsearchout'

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        else:               cmd = ['vsearch']
        # Other parameters
        cmd += ['--usearch_global', self.query,
                '-db',              self.db,
                '-blast6out',       self.out_path,
                '-threads',         self.cpus]
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return list(map(str, cmd))

    #-------------------------------- RUNNING --------------------------------#
    def run(self, verbose=True):
        """Simply run the VSEARCH search locally."""
        # Check the executable is available #
        if self.executable:
            self.executable.must_exist()
        else:
            from plumbing.check_cmd_found import check_cmd
            check_cmd('vsearch')
        # Create the output directory if it doesn't exist #
        self.out_path.directory.create_if_not_exists()
        # Optionally print the command #
        if verbose:
            print("Running VSEARCH command:\n    %s" % ' '.join(self.command))
        # Run it #
        cmd    = sh.Command(self.command[0])
        result = cmd(self.command[1:], _out=self._out, _err=self._err)
        # Return #
        return result

    #----------------------------- PARSE RESULTS -----------------------------#
    @property
    def results(self):
        """
        Parse the results and yield biopython SearchIO entries.

        Beware:
        Some databases are not unique on the id, and this causes the parser to
        complain about duplicate entries and raise exceptions such as:

            ValueError: The ID or alternative IDs of Hit 'DQ448783' exists
            in this QueryResult.

        Summary of the columns:
        https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

            qseqid sseqid pident length mismatch gapopen qstart qend sstart send
            evalue bitscore
        """
        with open(self.out_path, 'rt') as handle:
            for entry in SearchIO.parse(handle, 'blast-tab', ):
                yield entry

###############################################################################
class VSEARCHdb(FilePath):
    """A VSEARCH database one can search against."""

    def __repr__(self):
        return '<%s on "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, path):
        # Call parent constructor #
        FilePath.__init__(self, path)
        # Check if the corresponding FASTA exists #
        self.fasta_path = self.replace_extension('fasta')

    def __bool__(self):
        """Does the indexed database actually exist?"""
        return bool(self.replace_extension('udb'))

    def create_if_not_exists(self, *args, **kwargs):
        """If the indexed database has not been generated, generate it."""
        if not self: return self.makedb(*args, **kwargs)

    def makedb(self, output=None, stdout=None, verbose=False):
        # We need a corresponding fasta that exists #
        if not self.fasta_path:
            raise Exception("No fasta file at '%s'." % self.fasta_path)
        # Message #
        if verbose:
            msg = "Calling `makeudb_usearch` from '%s' to '%s'."
            print(msg % (self.fasta_path, self))
        # The output #
        if output is None: output = self
        # Options #
        options = ['--makeudb_usearch', self.fasta_path, '--output', output]
        # Call the program #
        sh.vsearch(*options, _out=stdout)
        # Return #
        return output
