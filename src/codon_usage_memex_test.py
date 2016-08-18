# Copyright (C) 2016 The President and Fellows of Harvard College
#
# This file may require additional software or modules to be installed to run
# run properly.
#
# The Genome Recoder software is available under an internal non-commercial
# research and academic use license.  Questions about this software or the
# licensing thereof can be addressed to  Office of Technology Development,
# Harvard University, email: otd@harvard.edu.
#
# @author Gleb Kuznetsov (kuznetsov@g.harvard.edu)

"""
Tests for codon_replacer.py.
"""

import unittest

from codon_usage_memex import build_codon_usage_dict
from codon_usage_memex import CodonUsageMemex

class TestCodonUsageMemex(unittest.TestCase):

    def test_build_codon_usage_dict(self):
        """Simple test to make sure the function at least runs without errors.
        """
        aa_to_codon_list_dict = build_codon_usage_dict(
                codon_usage_source_file='../data/ecoli-codon-usage.txt')


    def test_get_synonymous_codons(self):
        """Test getting synonymous codons.
        """
        # Re-assign random function.
        import random
        def fake_random():
            return 0.1
        random.random = fake_random

        CODON_USAGE_MEMEX = CodonUsageMemex(build_codon_usage_dict())

        print CODON_USAGE_MEMEX.get_synonymous_codons('TTA')



if __name__ == '__main__':
    unittest.main()

