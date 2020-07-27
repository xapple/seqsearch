# Built-in modules #

# Internal modules #
from autopaths.file_path import FilePath

# Third party modules #

# Constants #
tabular_keys = (
    'query_id',
    'subject_id',
    'perc_identity',
    'aln_length',
    'mismatch_count',
    'gap_open_count',
    'query_start',
    'query_end',
    'subject_start',
    'subject_end',
    'e_value',
    'bit_score',
)

###############################################################################
def parse_tabular(path):
    if not isinstance(path, FilePath): path = FilePath(path)
    for line in path:
        yield dict(zip(tabular_keys, line.split()))