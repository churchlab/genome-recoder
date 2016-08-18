"""
Tests for the homology module.
"""

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from biopython_util import InsertType
from homology import calc_homology
from homology import find_features_to_check_for_homology


class TestHomology(unittest.TestCase):
    """Tests for the homology module.
    """

    def test_calc_homology_same(self):
        """Test homology between identical sequences."""
        source_seq = 'ATTG'
        copy_seq = 'ATTG'
        homology_pair_obj = {
                'source_seq': source_seq,
                'copy_seq': copy_seq,
        }
        self.assertEqual(1.0, calc_homology(homology_pair_obj))


    def test_calc_homology_different(self):
        """Test homology between different sequences."""
        source_seq = 'ATTG'
        copy_seq = 'ATTC'
        homology_pair_obj = {
                'source_seq': source_seq,
                'copy_seq': copy_seq,
        }
        self.assertEqual(0.75, calc_homology(homology_pair_obj))


    def test_calc_different_lengths_raises_assertion_error(self):
        source_seq = 'ATTATATAT'
        copy_seq = 'ATTC'
        homology_pair_obj = {
                'source_seq': source_seq,
                'copy_seq': copy_seq,
        }
        self.assertRaises(AssertionError, calc_homology, homology_pair_obj)

class TestFindFeaturesToCheckForHomology(unittest.TestCase):

    def test_find_features_positive(self):
        before_junk = 'ATGCTAGCTAGCGGCGATAG'
        rbs = 'CCCATGTTT'
        overlap = 'ATGTTT'
        early_feature_1 = 'CTACTACTACTA'
        feature_1_seq = early_feature_1 + rbs + overlap
        rest_of_feature_2 = 'TTTCCCTTTCCC'
        feature_2_seq = overlap + rest_of_feature_2

        whole_seq = Seq(
                before_junk + feature_1_seq + rbs + overlap + rest_of_feature_2,
                generic_dna
        )
        seq_record = SeqRecord(whole_seq)

        ### Add all the features
        feature_1_start = len(before_junk)
        feature_1_end = len(before_junk) + len(feature_1_seq)
        feature_1_loc = FeatureLocation(feature_1_start, feature_1_end, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id='1')
        seq_record.features.append(feature_1)

        rbs_start = feature_1_end
        rbs_end = rbs_start + len(rbs)
        rbs_copy_loc = FeatureLocation(rbs_start, rbs_end, strand=1)
        rbs_copy_feature = SeqFeature(
                rbs_copy_loc, type=InsertType.FIX_OVERLAP_RBS_COPY, id='2')
        seq_record.features.append(rbs_copy_feature)

        head_copy_start = rbs_end
        head_copy_end = rbs_end + len(overlap)
        head_copy_loc = FeatureLocation(
                head_copy_start, head_copy_end, strand=1)
        head_copy_feature = SeqFeature(
                head_copy_loc, type=InsertType.FIX_OVERLAP_HEAD_COPY, id='3')
        seq_record.features.append(head_copy_feature)

        feature_2_start = rbs_end
        feature_2_end = rbs_end + len(feature_2_seq)
        feature_2_loc = FeatureLocation(feature_2_start, feature_2_end, strand=1)
        feature_2 = SeqFeature(feature_2_loc, type='CDS', id='4')
        seq_record.features.append(feature_2)


        ### Run the method in test.
        homologous_pairs = find_features_to_check_for_homology(seq_record)

        ### Handle asserts.
        self.assertEqual(2, len(homologous_pairs))


if __name__ == '__main__':
    unittest.main()
