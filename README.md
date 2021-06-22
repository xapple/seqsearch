[![PyPI version](https://badge.fury.io/py/seqsearch.svg)](https://badge.fury.io/py/seqsearch)

# `seqsearch` version 2.1.2

`seqsearch` is a python package for dealing sequence similarity searches (e.g. BLAST on DNA sequences) and automation.

It has many convenience methods that can automatically launch several types of search algorithms, as well as quick installation of sequence reference databases.

## Prerequisites

Since `seqsearch` is written in python, it is compatible with all operating systems: Linux, macOS and Windows. The only prerequisite is `python3` (which is often installed by default) along with the `pip3` package manager.

To check if you have `python3` installed, type the following on your terminal:

    $ python3 -V

If you do not have `python3` installed, please refer to the section [obtaining python3](docs/installing_tips.md#obtaining-python3).

To check you have `pip3` installed, type the following on your terminal:

    $ pip3 -V

If you do not have `pip3` installed, please refer to the section [obtaining pip3](docs/installing_tips.md#obtaining-pip3).

## Installing

To install the `seqsearch` package, simply type the following commands on your terminal:

    $ pip3 install --user seqsearch

Alternatively, if you want to install it for all users of the system:

    $ sudo pip3 install seqsearch

## Usage

Bellow are some examples to illustrate the various ways there are to use this package.

### Searches

You can parallelize BLAST searches by splitting the input into several files. It's easier to chop-up the input, because database chopping requires message passing across the nodes like mpiblast does (when and if it works).

Input chopping is fine as long as the database to search against fits in the RAM of the nodes. If the input is small and the database is large you can always switch them one for the other (in most cases).

    # This example is not completed yet. TODO.

### Databases

    # This example is not completed yet. TODO.

## Extra documentation

More documentation is available at:

<http://xapple.github.io/seqsearch/seqsearch>

This documentation is simply generated with:

    $ pdoc --html --output-dir docs --force seqsearch