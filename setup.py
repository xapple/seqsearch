#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Imports #
from setuptools import setup, find_namespace_packages
from os import path

# Load the contents of the README file #
this_dir = path.abspath(path.dirname(__file__))
readme_path = path.join(this_dir, 'README.md')
with open(readme_path, encoding='utf-8') as handle: readme = handle.read()

# Call setup #
setup(
    name             = 'seqsearch',
    version          = '2.0.6',
    description      = 'Sequence similarity searches (e.g. BLAST) made easy.',
    license          = 'MIT',
    url              = 'https://github.com/xapple/seqsearch',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_namespace_packages(exclude=['example']),
    install_requires = ['plumbing>=2.8.9', 'autopaths>=1.4.6', 'fasta>=2.1.0',
                        'biopython', 'sh', 'tqdm', 'ftputil', 'wget'],
    python_requires  = ">=3.6",
    long_description = readme,
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)