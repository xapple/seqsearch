#!/usr/bin/env python

"""
A script to test the "parallelblast" package.

Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this:
$ test.py query.fasta swissprot
"""

# The libraries we need #
import sys, os
from parallelblast import ParallelBLASTquery
# Get the shell arguments #
fa_path = sys.argv[1]
db_path = sys.argv[2]
# Check that the path is valid #
if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
# Params #
blast_params = {'algorithm': 'blastx' ,'-e': '1e-10', '-b': '0', '-v': '3', '-m': '7', '-a': 16}
slurm_params = {'time': '00:15:00', 'qos': 'short'}
split_size = '200K'
# Do it #
query = ParallelBLASTquery(fa_path, db_path, split_size, blast_params, slurm_params)
query.split()
query.run()
query.wait()
query.regroup()
# Check with non parallel version #
print "Diff between parallel and non parallel:" + query.check_np()