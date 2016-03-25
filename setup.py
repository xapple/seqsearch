from distutils.core import setup

setup(
      name             = 'seqsearch',
      version          = '1.0.1',
      description      = 'Sequence similarity searches made easy.',
      long_description = open('README.md').read(),
      license          = 'MIT',
      url              = 'https://github.com/xapple/seqsearch',
      author           = 'Lucas Sinclair',
      author_email     = 'lucas.sinclair@me.com',
      classifiers      = ['Topic :: Scientific/Engineering :: Bio-Informatics'],
      packages         = ['seqsearch'],
      requires         = ['plumbing', 'fasta', 'biopython', 'sh', 'shell_command', 'decorator', 'humanfriendly'],
)