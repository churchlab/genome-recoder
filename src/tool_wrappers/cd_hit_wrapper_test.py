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
