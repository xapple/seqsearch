#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import sys, os, multiprocessing, threading, shutil

# Internal modules #
from fasta import FASTA
from autopaths.tmp_path import new_temp_path
from autopaths.file_path import FilePath
from plumbing.cache import property_cached
from plumbing.slurm.job import JobSLURM

# Third party modules #
import sh
from ftputil import FTPHost
import Bio.Blast.NCBIXML

###############################################################################
class BLASTquery(object):
    """
    A blast job. Possibly the standard BLAST algorithm or
    BLASTP or BLASTX etc. Typically you could use it like this:

         import sys, os
         records_path = os.path.expanduser(sys.argv[1])
         centers_path = 'centers.fasta'
         db = parallelblast.BLASTdb(centers_path)
         db.makeblastdb()
         params = {'executable': "~/share/blastplus/blastn",
                   '-outfmt': 0,
                   '-evalue': 1e-2,
                   '-perc_identity': 97,
                   '-num_threads': 16}
         search = parallelblast.BLASTquery(records_path, db, params)
         search.run()

    You can also call search.non_block_run() to run maybe searches in parallel.
    """

    def __repr__(self): return '<%s object on %s>' % (self.__class__.__name__, self.query)

    def __init__(self, query_path, db_path,
                 seq_type     = 'prot' or 'nucl',     # The seq type of the query_path file
                 params       = None,                 # Add extra params for the command line
                 algorithm    = "blastn" or "blastp", # Will be autodetermined with seq_type
                 version      = "plus" or "legacy",   # Either blast+ or the old `blastall`
                 out_path     = None,                 # Where the results will be dropped
                 executable   = None,                 # If you want a specific binary give the path
                 cpus         = None,                 # The number of threads to use
                 num          = None,                 # When parallelizing, the number of this thread
                 slurm_params = None,                 # If you have special slurm parameters
                 _out         = None,                 # Store the stdout at this path
                 _err         = None):                # Store the stderr at this path
        # Save attributes #
        self.query        = FASTA(query_path)
        self.db           = BLASTdb(db_path, seq_type)
        self.seq_type     = seq_type
        self.version      = version
        self.algorithm    = algorithm
        self.num          = num
        self.params       = params if params else {}
        self.slurm_params = slurm_params if slurm_params else {}
        self._out         = _out
        self._err         = _err
        # Output #
        if out_path is None:         self.out_path = self.query.prefix_path + '.blastout'
        elif out_path.endswith('/'): self.out_path = out_path + self.query.prefix + '.blastout'
        else:                        self.out_path = out_path
        # Make it a file path #
        self.out_path = FilePath(self.out_path)
        # Executable #
        self.executable = FilePath(executable)
        # Cores to use #
        if cpus is None: self.cpus = multiprocessing.cpu_count()
        else:            self.cpus = cpus
        # Auto detect XML output #
        if self.out_path.extension == '.xml':
            if self.version == 'legacy': self.params['-xxxx']   = '5'
            if self.version == 'plus':   self.params['-outfmt'] = '5'

    @property
    def command(self):
        # Executable #
        if self.executable:            cmd = [self.executable.path]
        elif self.version == 'legacy': cmd = ["blastall", '-p', self.algorithm]
        else:                          cmd = [self.algorithm]
        # Essentials #
        if self.version == 'legacy':
            cmd += ['-d',  self.db, '-i',     self.query, '-o',   self.out_path, '-a',           self.cpus]
        if self.version == 'plus':
            cmd += ['-db', self.db, '-query', self.query, '-out', self.out_path, '-num_threads', self.cpus]
        # Options #
        for k,v in self.params.items(): cmd += [k, v]
        # Return #
        return map(str, cmd)

    #-------------------------------- RUNNING --------------------------------#
    def run(self):
        """Simply run the BLAST search locally."""
        out = self._out if self._out else '/dev/null'
        err = self._err if self._err else '/dev/null'
        sh.Command(self.command[0])(self.command[1:], _out=out, _err=err)
        if os.path.exists("error.log") and os.path.getsize("error.log") == 0: os.remove("error.log")

    def non_block_run(self):
        """Special method to run the query in a thread locally without blocking."""
        self.thread = threading.Thread(target=self.run)
        self.thread.daemon = True # So that they die when we die
        self.thread.start()

    def wait(self):
        """If you have run the query in a non-blocking way, call this method to pause
        until the query is finished."""
        try: self.thread.join(sys.maxint) # maxint timeout so that we can Ctrl-C them
        except KeyboardInterrupt: print("Stopped waiting on BLAST thread number %i" % self.num)

    @property_cached
    def slurm_job(self):
        """If you have access to a cluster with a SLURM queuing system, this property
        will return a SLURMjob object from which you can submit and check the status
        of the job."""
        return JobSLURM(' '.join(self.command), 'bash', self.out_path.directory, **self.slurm_params)

    #------------------------------- FILTERING -------------------------------#
    def filter(self, filtering):
        """We can do some special filtering on the results.
        For the moment only minimum coverage and minimum identity."""
        # Conditions #
        if 'min_coverage' in filtering and 'qcovs' not in self.params['-outfmt']:
            raise Exception("Can't filter on minimum coverage because it wasn't included.")
        if 'min_identity' in filtering and 'pident' not in self.params['-outfmt']:
            raise Exception("Can't filter on minimum identity because it wasn't included.")
        # Iterator #
        def filter_lines(blastout):
            cov_threshold = filtering.get('min_coverage', 0.0) * 100
            idy_threshold = filtering.get('min_identity', 0.0) * 100
            cov_position = self.params['-outfmt'].strip('"').split().index('qcovs') - 1
            idy_position = self.params['-outfmt'].strip('"').split().index('pident') - 1
            for line in blastout:
                coverage = float(line.split()[cov_position])
                identity = float(line.split()[idy_position])
                if coverage < cov_threshold: continue
                if identity < idy_threshold: continue
                else: yield line
        # Do it #
        temp_path = new_temp_path()
        with open(temp_path, 'w') as handle: handle.writelines(filter_lines(self.out_path))
        os.remove(self.out_path)
        shutil.move(temp_path, self.out_path)

    @property
    def results(self):
        """Parse the results."""
        # Check for XML #
        if self.params.get('-outfmt', 0) == '5':
            with open(self.out_path, 'rb') as handle:
                for entry in Bio.Blast.NCBIXML.parse(handle): yield entry
        # Default case #
        else:
            for line in self.out_path: yield line.split()

###############################################################################
class BLASTdb(FASTA):
    """A BLAST database one can search against."""

    def __repr__(self): return '<%s on "%s">' % (self.__class__.__name__, self.path)

    def __init__(self, fasta_path, seq_type='nucl' or 'prot'):
        if hasattr(fasta_path, 'seq_type'): self.seq_type = fasta_path.seq_type
        else:                               self.seq_type = seq_type
        FASTA.__init__(self, fasta_path)

    def makeblastdb(self, logfile=None, out=None):
        # Message #
        print("Calling `makeblastdb` on '%s'..." % self)
        # Options #
        options = ['-in', self.path, '-dbtype', self.seq_type]
        # Add a log file #
        if logfile is not None: options += ['-logfile', logfile]
        # Call the program #
        if out is not None: sh.makeblastdb(*options, _out=str(out))
        else:               sh.makeblastdb(*options)

###############################################################################
def install_blast(base_dir):
    """Deprecated, look into 'home_linux/setup/bioinfo_tools'
    On Ubuntu: sudo apt-get install ncbi-blast+
    """
    # Default location #
    if base_dir is None: base_dir = os.environ.get('HOME', '/') + '/programs/blast/'
    # Download from FTP #
    ftp_url = "ftp.ncbi.nlm.nih.gov"
    ftp_dir = "/blast/executables/blast+/LATEST/"
    pattern = 'ncbi-blast-*+-src.zip'
    ftp = FTPHost(ftp_url, "anonymous")
    ftp.chdir(ftp_dir)
    files = ftp.listdir(self.ftp.curdir)
    ftp.download(source, dest)