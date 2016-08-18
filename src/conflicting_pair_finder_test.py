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
Tests for methods that find conflicting pairs (overlaps and rbs conflicts).
"""

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from conflicting_pair_finder import is_potential_rbs_conflict


class TestFindConflictingFeaturePairs(unittest.TestCase):
    """Test methods that find feature pairs that are either overlapping
    or close enough to potentially affect each other's RBS.
    """

    def test_is_potential_rbs_conflict(self):
        feature_1_seq = 'ATGTTTAAACCC'
        gap = 'CCC'
        feature_2_seq = 'ATGTTTAAACCC'
        whole_seq = feature_1_seq + gap + feature_2_seq

        # Same direction.
        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
        self.assertTrue(is_potential_rbs_conflict(feature_1, feature_2))

        # Reverse same direction.
        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
        self.assertTrue(is_potential_rbs_conflict(feature_1, feature_2))

        # Opposite dirs.
        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
        self.assertTrue(is_potential_rbs_conflict(feature_1, feature_2))

        # Facing each other, false.
        # NOTE: Safe to assume actual features will be bigger than what they
        # are in this test.
        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=-1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')
        self.assertFalse(is_potential_rbs_conflict(feature_1, feature_2))

        # Not both CDS, false.
        # NOTE: Safe to assume actual features will be bigger than what they
        # are in this test.
        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='misc_RNA', id='2')
        self.assertFalse(is_potential_rbs_conflict(feature_1, feature_2))



    def test_is_potential_rbs_conflict__false_too_far(self):
        feature_1_seq = 'ATGTTTAAACCC'
        too_large_gap = 'AAATTTCCCAAATTTCCCAAATTTCCCAAATTTCCCAAATTTCCC'
        feature_2_seq = 'ATGTTTAAACCC'
        whole_seq = feature_1_seq + too_large_gap + feature_2_seq

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')

        feature_2_loc = FeatureLocation(
                len(whole_seq) - len(feature_2_seq),
                len(whole_seq),
                strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='2')

        self.assertFalse(is_potential_rbs_conflict(feature_1, feature_2))


if __name__ == '__main__':
    unittest.main()
