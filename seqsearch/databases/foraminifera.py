# Built-in modules #
import os, tarfile

# First party modules #
from seqsearch.databases import Database
from fasta import FASTA
from plumbing.autopaths import AutoPaths, FilePath

# Third party modules #

# Constants #
home = os.environ.get('HOME', '~') + '/'

###############################################################################
class Foraminifera(Database):
    """This is a custom database containing exlcusively Foraminifera sequences.

    https://genev.unige.ch/research/laboratory/Jan-Pawlowski

    You should place the file "foram_db_flo.fasta" in:  ~/databases/foraminifera/
    Then you can run this:
    
            from seqsearch.databases.foraminifera import foraminifera

    """

    short_name = "foraminifera"

    all_paths = """
    /foram_db_flo.fasta
    /foram_mothur.fasta
    /foram_mothur.tax
    """

    def __init__(self, base_dir=None):
        # Base directory #
        if base_dir is None: base_dir = home
        self.base_dir = base_dir + 'databases/' + self.short_name + '/'
        self.p        = AutoPaths(self.base_dir, self.all_paths)
        # The results #
        self.alignment = FASTA(self.p.mothur_fasta)
        self.taxonomy  = FilePath(self.p.mothur_tax)

    def process(self):
        # The file that was received by email T_T #
        raw = FASTA(self.p.flo)
        # Open files #
        self.alignment.create()
        self.taxonomy.create()
        # Loop #
        for seq in raw:
            # Parse #
            name = seq.id[11:].split('|')
            num  = name.pop(0)
            # Check #
            for x in name: assert ';' not in x
            for x in name: assert '\t' not in x
            # Make ranks #
            ranks = ['Eukaryota'                       , # Domain
                     'Foraminifera'                    , # Phylum
                     name[0]                           , # Class
                     name[1]                           , # Order
                     name[2] + " (" + name[3] + ")"    , # Family
                     name[4]                           , # Genus
                     name[5]]                            # Species
            # The taxonomy string #
            tax_line = ';'.join(ranks)
            # Add sequence to the new fasta file #
            self.alignment.add_str(str(seq.seq), name=num)
            # Add the taxonomy to the tax file #
            self.taxonomy.add_str(num + '\t' + tax_line + '\n')
        # Close files #
        self.alignment.close()
        self.taxonomy.close()

###############################################################################
foraminifera = Foraminifera()