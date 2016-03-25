[![changelog](http://allmychanges.com/p/python/seqsearch/badge/)](http://allmychanges.com/p/python/seqsearch/?utm_source=badge) [![PyPI version](https://badge.fury.io/py/seqsearch.svg)](https://badge.fury.io/py/seqsearch)

# `seqsearch` version 1.0.3

Sequence similarity searches made easy

### Introduction
You can parallelize BLAST searches by splitting the input and submitting jobs to a SLURM cluster.

The blast jobs are run on a SLURM queue with parallelization by input chopping. Not database chopping which requires message passing across the nodes like mpiblast does (when and if it works).
Input chopping is fine as long as the database to search against fits in the RAM of the nodes. If the input is small and the database is large you can always switch them one for the other (in most cases) and apply input splitting.

For a usage example look at the `test.py` file.

### Code documentation
More documentation is available at:
http://xapple.github.io/seqsearch/