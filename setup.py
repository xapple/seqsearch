from distutils.core import setup

setup(
    name             = 'seqsearch',
    version          = '1.0.3',
    description      = 'Sequence similarity searches made easy.',
    license          = 'MIT',
    url              = 'https://github.com/xapple/seqsearch',
    author           = 'Lucas Sinclair',
    author_email     = 'lucas.sinclair@me.com',
    classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
    packages         = ['seqsearch'],
    install_requires = ['plumbing', 'fasta', 'biopython', 'sh', 'shell_command', 'decorator', 'humanfriendly'],
    long_description = open('README.md').read(),
)