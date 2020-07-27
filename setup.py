from setuptools import setup, find_packages

setup(
    name             = 'seqsearch',
    version          = '1.2.2',
    description      = 'Sequence similarity searches (e.g. BLAST) made easy.',
    license          = 'MIT',
    url              = 'https://github.com/xapple/seqsearch',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = find_packages(),
    install_requires = ['plumbing>=2.8.5', 'autopaths>=1.4.4', 'fasta>=2.0.9',
                        'biopython', 'decorator', 'sh', 'tqdm', 'ftputil',
                        'wget'],
    long_description = open('README.md').read(),
    long_description_content_type = 'text/markdown',
    include_package_data = True,
)