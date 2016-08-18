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
Tests for the CD-Hit wrapper module.
"""

import unittest

from cd_hit_wrapper import CDHitWrapper


class TestCDHitWrapper(unittest.TestCase):
    """Tests for the wrapper object."""

    def test_cd_hit_wrapper(self):
        WRAPPER_OBJ = CDHitWrapper()
        SEQUENCES = [
                'AAAAAAAAAAAAAAA',
                'TTTTTTTTTTTTTTT',
                'CCCCCCCCCCCCCCC',
                'AAAAAAAAAAAAAAAAAA'
        ]
        high_homology_sets = (
                WRAPPER_OBJ.detect_high_homology_pairs(SEQUENCES))

        EXPECTED_HIGH_HOMOLOGY_SETS = [set([0,3])]

        # TODO: Make this pass.
        self.assertItemsEqual(EXPECTED_HIGH_HOMOLOGY_SETS, high_homology_sets)


if __name__ == '__main__':
    unittest.main()
