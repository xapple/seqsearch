#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Written by Lucas Sinclair.
MIT Licensed.
Contact at www.sinclair.bio
"""

# Built-in modules #
import sh

# First party modules #
from seqsearch.databases import Database
from autopaths.auto_paths import AutoPaths
from autopaths.dir_path   import DirectoryPath
from plumbing.cache import property_cached
from fasta import FASTA

###############################################################################
class Pfam(Database):
    """The Pfam database is a large collection of protein families,
    each represented by multiple sequence alignments and HMMs.

    http://pfam.xfam.org

    Pfam-A is the manually curated portion of the database that
    contains over 16,000 entries. An automatically generated supplement
    is provided called Pfam-B. Pfam-B was discontinued.

    To install:
        from seqsearch.databases.pfam import pfam
        pfam.download()
        pfam.ungzip()

    It will put it in ~/databases/pfam
    """

    all_paths = """
    /raw/
    /unzipped/Pfam-A.hmm
    /unzipped/Pfam-A.fasta
    /unzipped/Pfam-A.seed
    /specific/
    """

    short_name = "pfam"
    ftp_url    = "ftp.ebi.ac.uk"
    ftp_dir    = "/pub/databases/Pfam/current_release/"
    files      = ("Pfam-A.hmm.gz", "Pfam-A.fasta.gz", "Pfam-A.seed.gz")

    @property_cached
    def hmm_db(self):
        hmm_db = self.autopaths.hmm
        hmm_db.seq_type = 'hmm_prot'
        return hmm_db

    @property_cached
    def fasta(self):
        fasta = FASTA(self.autopaths.fasta)
        return fasta

    @property_cached
    def seeds(self):
        seeds = FASTA(self.autopaths.seed)
        return seeds

###############################################################################
pfam = Pfam("hmm_prot")

###############################################################################
class SpecificFamily(object):
    """When you are interested in having an HMM 'database' with only
    one specific Pfam in it. As well as the associated proteins that
    are part of it."""

    all_paths = """
    /model.hmm
    /proteins.fasta
    /subsampled.fasta
    """

    def __init__(self, fam_name):
        self.fam_name = fam_name
        self.base_dir = DirectoryPath(pfam.autopaths.specific_dir + self.fam_name)
        self.p        = AutoPaths(self.base_dir, self.all_paths)

    @property_cached
    def hmm_db(self):
        """Create an HMM 'database' with only one profile in it."""
        hmm_db = self.p.model
        hmm_db.seqtype = 'hmm_prot'
        if not hmm_db.exists:
            print(sh.hmmfetch('-o', hmm_db, pfam.hmm_db, self.fam_name))
            assert hmm_db
        return hmm_db

    @property_cached
    def fasta(self):
        """Make a fasta file with all uniprot proteins that are related to
        this family."""
        fasta = FASTA(self.p.proteins)
        if not fasta.exists:
            fasta.create()
            for seq in pfam.fasta:
                if self.fam_name in seq.description: fasta.add_seq(seq)
            fasta.close()
            assert fasta
        # Return #
        return fasta

    @property_cached
    def subsampled(self):
        subsampled = FASTA(self.p.subsampled)
        if not subsampled.exists:
            self.fasta.subsample(down_to=30, new_path=subsampled)
            self.add_taxonomy(subsampled)
        return subsampled

    def add_taxonomy(self, fasta):
        """Add taxonomic information to the fastas file"""
        ids        = [seq.description for seq in fasta]
        accesions  = [seq.description.split()[1] for seq in fasta]
        taxonomies = map(self.uniprot_acc_to_taxonmy, accesions)
        naming_dict  = {ids[i]: ids[i] + taxonomies[i] for i in range(len(ids))}
        fasta.rename_sequences(naming_dict, in_place=True)

    def uniprot_acc_to_taxonmy(self, accesion):
        """From one uniprot ID to taxonomy"""
        from bioservices import UniProt
        u = UniProt()
        data = u.search(accesion, frmt="xml")
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(data, "html.parser")
        return ' (' + ', '.join([t.text for t in soup.find_all('taxon')]) + ')'
