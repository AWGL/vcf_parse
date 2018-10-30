#!/anaconda3/envs/python2/bin/python

"""
known_variants.py

Object that deals with loading and annotating known variants. 
Loaded as part of the vcf_parse.py program.

Author:     Erik Waskiewicz
Created:    30 Oct 2018
Version:    0.1.0
Updated:    30 Oct 2018
"""


import logging


# -- KNOWN VARIANTS CLASS ---------------------------------------------

class known_variants:
    def __init__(self):
        """
        Object properties that are loaded when the oject is created.
        The logger deals with all status messages from the object and
        is a child of the main logger vcf_parse.
        """
        self.logger = logging.getLogger('vcf_parse.known')


    def load_known_variants(self, in_vcf):
        """
        check that there is a Classification column
        Load in vcf and save as list
        """
        pass


    def apply_known_variants(self, bedfile, in_vcf, out_folder):
        """
        Compares variant id with lost of known variants, annotates if there is a match
        """
        pass
 