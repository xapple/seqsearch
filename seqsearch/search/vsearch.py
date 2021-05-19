#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #

# First party modules #

# Internal modules #
from fasta import FASTA
from seqsearch.search.core import CoreSearch

# Third party modules #
import sh

###############################################################################
class VSEARCHquery(CoreSearch):
    """
    A vsearch job.

    Example commands:

        vsearch --usearch_global queries.fas --db queries.fas --id 0.6
                --blast6out results.blast6 --maxaccepts 0 --maxrejects 0

        ./vsearch --usearch_global queries.fsa --db database.fsa --id 0.9
                  --alnout alnout.txt

    The interesting options form the manual:

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

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [self.executable.path]
        else:               cmd = ['vsearch']
        # Other parameters
        cmd += ['--usearch_global', self.query,
                '-db',              self.db,
                '-maxaccepts',      20,
                '-blast6out',       self.out_path,
                '-threads',         self.cpus]
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return list(map(str, cmd))

    def run(self, verbose=False):
        """Simply run the VSEARCH search locally."""
        # Optionally print the command #
        if verbose:
            print("Running VSEARCH command:\n    %s" % ' '.join(self.command))
        # Run it #
        sh.Command(self.command[0])(self.command[1:],
                                    _out=self._out,
                                    _err=self._err)

###############################################################################
class VSEARCHdb(FASTA):
    """A VSEARCH database one can search against."""

    def __repr__(self):
        return '<%s on "%s">' % (self.__class__.__name__, self.path)

    def __bool__(self):
        """Does the indexed database actually exist?"""
        return bool(self.replace_extension('udb'))

    def create_if_not_exists(self, *args, **kwargs):
        """If the indexed database has not been generated, generate it."""
        if not self: return self.makedb(*args, **kwargs)

    def makedb(self, output=None, stdout=None, verbose=False):
        # Message #
        if verbose: print("Calling `makeudb_usearch` on '%s'..." % self)
        # The output #
        if output is None: output = self.replace_extension('udb')
        # Options #
        options = ['--makeudb_usearch', self.path, '--output', output]
        # Call the program #
        sh.vsearch(*options, _out=stdout)
