#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import os, shutil

# First party modules #
from fasta import FASTA
from autopaths.tmp_path  import new_temp_path

# Internal modules #
from seqsearch.search.core import CoreSearch

# Third party modules #
import Bio.Blast.NCBIXML
from Bio import SearchIO

# Module for launching shell commands #
if os.environ.get('INSIDE_PYCHARM'): import pbs3 as sh
else:                                import sh

###############################################################################
class BLASTquery(CoreSearch):
    """
    A blast job. Possibly the standard BLAST algorithm or
    BLASTP or BLASTX etc. Typically you could use it like this:

         import sys, os
         records_path = os.path.expanduser(sys.argv[1])
         centers_path = 'centers.fasta'
         db = seqsearch.BLASTdb(centers_path)
         db.makeblastdb()
         params = {'executable':    "~/share/blastplus/blastn",
                   '-outfmt':        0,
                   '-evalue':        1e-2,
                   '-perc_identity': 97,
                   '-num_threads':   16}
         search = seqsearch.BLASTquery(records_path, db, params)
         search.run()

    You can also call search.non_block_run() to run many searches in parallel.
    """

    extension = 'blastout'

    def __init__(self, *args, **kwargs):
        # Parent constructor #
        super(BLASTquery, self).__init__(*args, **kwargs)
        # Auto detect XML output #
        if self.out_path.extension == '.xml': self.params['-outfmt'] = '5'
        # The database to search against #
        self.db = BLASTdb(self.db, self.seq_type)

    @property
    def command(self):
        # Executable #
        if self.executable: cmd = [str(self.executable)]
        else:               cmd = [self.algorithm]
        # Other parameters
        cmd += ['-db',          self.db,
                '-query',       self.query,
                '-out',         self.out_path,
                '-num_threads', self.cpus]
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return list(map(str, cmd))

    #-------------------------------- RUNNING --------------------------------#
    def run(self, verbose=False):
        """Simply run the BLAST search locally."""
        # Check the executable is available #
        if self.executable:
            self.executable.must_exist()
        else:
            from plumbing.check_cmd_found import check_cmd
            check_cmd(self.command[0])
        # Create the output directory if it doesn't exist #
        self.out_path.directory.create_if_not_exists()
        # Optionally print the command #
        if verbose:
            print("Running BLAST command:\n    %s" % ' '.join(self.command))
        # Run it #
        cmd    = sh.Command(self.command[0])
        result = cmd(self.command[1:], _out=self._out, _err=self._err)
        # Clean up #
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0:
            os.remove("error.log")
        # Return #
        return result

    #------------------------------- FILTERING -------------------------------#
    def filter(self, filtering):
        """
        We can do some special filtering on the results.
        For the moment only minimum coverage and minimum identity.
        """
        # Conditions #
        if 'min_coverage' in filtering and \
           'qcovs' not in self.params['-outfmt']:
            msg = "Can't filter on minimum coverage because it wasn't included."
            raise Exception(msg)
        if 'min_identity' in filtering and \
           'pident' not in self.params['-outfmt']:
            msg = "Can't filter on minimum identity because it wasn't included."
            raise Exception(msg)
        # Iterator #
        def filter_lines(blastout):
            cov_threshold = filtering.get('min_coverage', 0.0) * 100
            idy_threshold = filtering.get('min_identity', 0.0) * 100
            outfmt_str    = self.params['-outfmt'].strip('"').split()
            cov_position  = outfmt_str.index('qcovs')  - 1
            idy_position  = outfmt_str.index('pident') - 1
            for line in blastout:
                coverage = float(line.split()[cov_position])
                identity = float(line.split()[idy_position])
                if coverage < cov_threshold: continue
                if identity < idy_threshold: continue
                else:yield line
        # Do it #
        temp_path = new_temp_path()
        with open(temp_path, 'w') as handle:
            handle.writelines(filter_lines(self.out_path))
        os.remove(self.out_path)
        shutil.move(temp_path, self.out_path)

    #----------------------------- PARSE RESULTS -----------------------------#
    @property
    def results(self):
        """Parse the results and yield biopython SearchIO entries."""
        # Get the first number of the outfmt #
        outfmt_str = self.params.get('-outfmt', '0').strip('"').split()
        number = outfmt_str[0]
        # Check for XML #
        if number == '5':
            with open(self.out_path, 'rb') as handle:
                for entry in Bio.Blast.NCBIXML.parse(handle):
                    yield entry
        # Check for tabular #
        elif number == '6':
            with open(self.out_path, 'rt') as handle:
                for entry in SearchIO.parse(handle, 'blast-tab'):
                    yield entry
        # Check for tabular with comments #
        elif number == '7':
            with open(self.out_path, 'rt') as handle:
                for entry in SearchIO.parse(handle, 'blast-tab', comments=True):
                    yield entry
        # Default case #
        else:
            for line in self.out_path:
                yield line.split()

###############################################################################
class BLASTdb(FASTA):
    """A BLAST database one can search against."""

    def __repr__(self):
        return '<%s on "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, fasta_path, seq_type='nucl' or 'prot'):
        # Check if the FASTA already has the seq_type set #
        if hasattr(fasta_path, 'seq_type'): self.seq_type = fasta_path.seq_type
        else:                               self.seq_type = seq_type
        # Call parent constructor #
        FASTA.__init__(self, fasta_path)

    def __bool__(self):
        """Does the indexed database actually exist?"""
        return bool(self + '.nsq')

    def create_if_not_exists(self, *args, **kwargs):
        """If the indexed database has not been generated, generate it."""
        if not self: return self.makedb(*args, **kwargs)

    def makedb(self, logfile=None, stdout=None, verbose=False):
        # Message #
        if verbose: print("Calling `makeblastdb` on '%s'..." % self)
        # Options #
        options = ['-in', self.path, '-dbtype', self.seq_type]
        # Add a log file #
        if logfile is not None: options += ['-logfile', logfile]
        # Call the program #
        sh.makeblastdb(*options, _out=stdout)
