# seqsearch version 1.0.0

Sequence similarity searches made easy

### Introduction
You can parallelize BLAST searches by splitting the input and submitting jobs to a SLURM cluster.

The blast jobs are run on a SLURM queue with parallelization by input chopping. Not database chopping which requires message passing across the nodes like mpiblast does (when and if it works).
Input chopping is fine as long as the database to search against fits in the RAM of the nodes. If the input is small and the database is large you can always switch them one for the other (in most cases) and apply input splitting.

For a usage example look at the test.py file.