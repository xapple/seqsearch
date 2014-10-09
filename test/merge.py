"""
Copied and adapted from https://bitbucket.org/peterjc/galaxy-central/src/5cefd5d5536e/tools/ncbi_blast_plus/blast.py

Tested working with BLAST 2.2.28+

More info at: http://www.biostars.org/p/53946/
"""

###############################################################################
def merge_blast_xml(inputs, output):
    """Will merge the results from several BLAST searches
    when XML output is used."""
    for i,handle in enumerate(inputs):
        path = handle.name
        header = handle.readline()
        if not header:
            raise Exception("BLAST XML file '%s' was empty" % path)
        if header.strip() != '<?xml version="1.0"?>':
            raise Exception("BLAST file '%s' is not an XML file" % path)
        line = handle.readline()
        header += line
        if line.strip()[0:59] != '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"':
            raise Exception("BLAST file '%s' is not a valid XML file" % path)
        while True:
            line = handle.readline()
            if not line:
                raise Exception("BLAST XML file '%s' ended prematurely" % path)
            header += line
            if "<Iteration>" in line: break
            if len(header) > 10000:
                raise Exception("BLAST file '%s' has a too long a header" % path)
        if "<BlastOutput>" not in header:
            raise Exception("BLAST XML file '%s' header's seems bad" % path)
        if i == 0:
            output.write(header)
            old_header = header
        elif old_header[:300] != header[:300]:
            raise Exception("BLAST XML headers don't match" % path)
        else: output.write("    <Iteration>\n")
        for line in handle:
            if "</BlastOutput_iterations>" in line: break
            output.write(line)
    output.write("  </BlastOutput_iterations>\n")
    output.write("</BlastOutput>\n\n")
    output.flush()