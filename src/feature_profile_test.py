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
Tests for feature_profile.py
"""

import unittest

from Bio.Alphabet import generic_dna
from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from codon_usage_memex import CodonUsageMemex
from feature_profile import FeatureProfile
from feature_profile import CodonRarityFeatureProfile
from feature_profile import GCContentFeatureProfile
from feature_profile import SecondaryStructureFeatureProfile


class TestFeatureProfile(unittest.TestCase):
    """Tests for the default base class.
    """

    def setUp(self):
        """Override."""
        # Okay to override in this test, but in real code you should
        # actually define the attribute in any sub-class.
        FeatureProfile._name = 'base'


    def test_base_profile(self):
        """Test the base class implementation.
        """
        feature_1_seq = 'ATGTTTGGG'
        other = 'TTTAAACCCTTT'
        whole_seq = feature_1_seq + other
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(0, len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Test profile construction.
        profile = FeatureProfile(feature_1.id, seq_record)
        self.assertEqual(len(feature_1_seq) / 3, len(profile.values))

        # Test individual member methods.
        self.assertEqual(0, profile.codon_index_to_seq_record_pos(0))
        self.assertEqual(3, profile.codon_index_to_seq_record_pos(1))
        self.assertEqual(6, profile.codon_index_to_seq_record_pos(2))


    def test_constructor__forward(self):
        """Make sure that we get the polarity_aware_* attributes correct."""
        before = 'CCC'
        feature_1_seq = 'ATGTTTGGG'
        after = 'TTTAAACCCTTT'
        whole_seq = before + feature_1_seq + after
        SEQ = Seq(whole_seq, generic_dna)
        SEQ_RECORD_POSITIVE = SeqRecord(SEQ)

        feature_1_loc = FeatureLocation(
                len(before), len(before) + len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        SEQ_RECORD_POSITIVE.features.append(feature_1)
        profile = FeatureProfile(feature_1.id, SEQ_RECORD_POSITIVE)
        underlying_feature_seq = str(profile.polarity_aware_underlying_seq[
                profile.polarity_aware_feature_location_start:
                profile.polarity_aware_feature_location_end])
        self.assertEqual(feature_1_seq, underlying_feature_seq)


    def test_constructor__negative(self):
        """Make sure that we get the polarity_aware_* attributes correct."""
        before = 'CCC'
        feature_1_seq = 'ATGTTTGGG'
        after = 'TTTAAACCCTTT'
        whole_seq = before + reverse_complement(feature_1_seq) + after
        SEQ = Seq(whole_seq, generic_dna)
        SEQ_RECORD_NEGATIVE = SeqRecord(SEQ)

        feature_1_loc = FeatureLocation(
                len(before), len(before) + len(feature_1_seq), strand=-1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        SEQ_RECORD_NEGATIVE.features.append(feature_1)
        profile = FeatureProfile(feature_1.id, SEQ_RECORD_NEGATIVE)

        EXPECTED_FEATURE_START = len(after)
        EXPECTED_FEATURE_END = EXPECTED_FEATURE_START + len(feature_1_seq)
        self.assertEqual(EXPECTED_FEATURE_START,
                profile.polarity_aware_feature_location_start)
        self.assertEqual(EXPECTED_FEATURE_END,
                profile.polarity_aware_feature_location_end)

        underlying_feature_seq = str(profile.polarity_aware_underlying_seq[
                profile.polarity_aware_feature_location_start:
                profile.polarity_aware_feature_location_end])
        self.assertEqual(feature_1_seq, underlying_feature_seq)


    def test_get_window_seq__positive(self):
        """Tests getting the window sequence for computing a profile value
        on the positive strand.
        """
        # Build the sequence.
        BEFORE = 'CCCTTTGGG'
        FEATURE_1_SEQ = 'ATGTTTGGG'
        AFTER = 'TTTAAACCCTTT'
        WHOLE_SEQ = BEFORE + FEATURE_1_SEQ + AFTER
        SEQ = Seq(WHOLE_SEQ, generic_dna)
        SEQ_RECORD = SeqRecord(SEQ)

        # Build and add the feature.
        FEATURE_1_LOC = FeatureLocation(
                len(BEFORE), len(BEFORE) + len(FEATURE_1_SEQ), strand=1)
        FEATURE_1 = SeqFeature(FEATURE_1_LOC, type='CDS', id=1)
        SEQ_RECORD.features.append(FEATURE_1)

        # Build the profile.
        SLIDING_WINDOW_SIZE = 6
        PROFILE_KWARGS = {
                'sliding_window_size': SLIDING_WINDOW_SIZE
        }
        PROFILE = FeatureProfile(FEATURE_1.id, SEQ_RECORD, **PROFILE_KWARGS)

        EXPECTED_WINDOW_SEQ = BEFORE[-3:] + FEATURE_1_SEQ[:3]
        self.assertEqual(EXPECTED_WINDOW_SEQ,
                PROFILE.get_window_seq(0, FEATURE_1_SEQ)
        )


    def test_get_window_seq__negative(self):
        """Tests getting the window sequence for computing a profile value
        on the negative strand.
        """
        # Build the sequence.
        BEFORE = 'CCCTTTGGG'
        FEATURE_1_SEQ = 'ATGTTTGGG'
        AFTER = 'ATGCGATTGGATCGA'
        WHOLE_SEQ = BEFORE + reverse_complement(FEATURE_1_SEQ) + AFTER
        SEQ = Seq(WHOLE_SEQ, generic_dna)
        SEQ_RECORD = SeqRecord(SEQ)

        # Build and add the feature.
        FEATURE_1_LOC = FeatureLocation(
                len(BEFORE), len(BEFORE) + len(FEATURE_1_SEQ), strand=-1)
        FEATURE_1 = SeqFeature(FEATURE_1_LOC, type='CDS', id=1)
        SEQ_RECORD.features.append(FEATURE_1)

        # Build the profile.
        SLIDING_WINDOW_SIZE = 6
        PROFILE_KWARGS = {
                'sliding_window_size': SLIDING_WINDOW_SIZE
        }
        PROFILE = FeatureProfile(FEATURE_1.id, SEQ_RECORD, **PROFILE_KWARGS)

        EXPECTED_WINDOW_SEQ = reverse_complement(WHOLE_SEQ)[
                len(AFTER) - 3: len(AFTER) + 3]
        self.assertEqual(EXPECTED_WINDOW_SEQ,
                PROFILE.get_window_seq(0, FEATURE_1_SEQ)
        )


class TestGCContentFeatureProfile(unittest.TestCase):
    """Tests for the GC content profiler.
    """

    def test_gc_profile(self):
        """Basic test.
        """
        # Override sliding_window_size for testing.
        GCContentFeatureProfile.sliding_window_size = 6

        # Create the test sequence.
        before_feature = 'CCC'
        feature_1_seq = 'AAATTTGGGCCC'
        after_feature = 'GGG'
        whole_seq = before_feature + feature_1_seq + after_feature
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        # Create the feature object.
        feature_1_start = len(before_feature)
        feature_1_end = len(before_feature) + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start, feature_1_end, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Create the profile.
        profile = GCContentFeatureProfile(feature_1.id, seq_record)

        EXPECTED_PROFILE_SCORES = [0.5, 0, 0.5, 1.0]

        self.assertEqual(len(EXPECTED_PROFILE_SCORES), len(profile.values))
        self.assertEqual(EXPECTED_PROFILE_SCORES, profile.values)


    def test_reverse_strand(self):
        """Make sure that the reverse strand score is actually calculated
        in the correct direction.
        """
        # Override sliding_window_size for testing.
        GCContentFeatureProfile.sliding_window_size = 6

        # Create the test sequence.
        before_feature = 'CCCGGG'
        feature_seq = 'GGGCCCGGGAAA'
        after_feature = 'AAATTT'
        whole_seq = before_feature + feature_seq + after_feature
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        # Create the forward feature object.
        feature_start = len(before_feature)
        feature_end = len(before_feature) + len(feature_seq)
        forward_feature_loc = FeatureLocation(
                feature_start, feature_end, strand=1)
        forward_feature = SeqFeature(forward_feature_loc, type='CDS', id=1)
        seq_record.features.append(forward_feature)

        # Create the reverse feature.
        reverse_feature_loc = FeatureLocation(
                feature_start, feature_end, strand=-1)
        reverse_feature = SeqFeature(reverse_feature_loc, type='CDS', id=2)
        seq_record.features.append(reverse_feature)

        # Create the profiles.
        forward_feature_profile = GCContentFeatureProfile(
                forward_feature.id, seq_record)
        reverse_feature_profile = GCContentFeatureProfile(
                reverse_feature.id, seq_record)

        EXPECTED_FORWARD_PROFILE = [1.0, 1.0, 1.0, 0.5]
        EXPECTED_REVERSE_PROFILE = [0.0, 0.5, 1.0, 1.0]

        self.assertEqual(EXPECTED_FORWARD_PROFILE,
                forward_feature_profile.values)
        self.assertEqual(EXPECTED_REVERSE_PROFILE,
                reverse_feature_profile.values)


class TestSecondaryStructureFeatureProfile(unittest.TestCase):
    """Tests for the GC content profiler.
    """

    def test_secondary_structure_profile(self):
        """Basic test.
        """
        # Create the test sequence.
        before_feature = 'GAGCTAGTGCGCTGATAGCCGGATAGCGCAAAATTGCTAGATAGCGATCGATCGAGAGAGAGATCTGAAATATATTTAGAG'
        feature_1_seq = 'ATGTTTGGGCCCGATGCTAATCCGAATCGACTA'
        after_feature = 'GGGATGTTTCGGGAGGGGTAGGATAGATAGAGGATATATAG'
        whole_seq = before_feature + feature_1_seq + after_feature
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        # Create the feature object.
        feature_1_start = len(before_feature)
        feature_1_end = len(before_feature) + len(feature_1_seq)
        feature_1_loc = FeatureLocation(
                feature_1_start, feature_1_end, strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Create the profile.
        profile = SecondaryStructureFeatureProfile(feature_1.id, seq_record)

        # NOTE: Expected sequence pre-caculated. The point is just to make
        # sure nothing inadvertendly changes in the algorithm. If a change is
        # necessary, please update the test.
        EXPECTED_PROFILE_VALUES = [
                -8.5, -7.8, -6.9, -5.6, -5.2, -6.4, -7.4, -5.2, -4.0, -4.5, -4.5
        ]
        self.assertEqual(len(EXPECTED_PROFILE_VALUES), len(profile.values))
        self.assertEqual(EXPECTED_PROFILE_VALUES, profile.values)


class TestCodonRarityFeatureProfile(unittest.TestCase):
    """Tests for the CodonRarityFeatureProfile.
    """

    def test_codon_rarity_profile(self):
        # Override sliding_window_size for testing.
        CodonRarityFeatureProfile.sliding_window_size = 15

        before = 'TTTAAACCCTTTGGG'
        feature_1_seq = 'ATGTTTGGG'
        whole_seq = before + feature_1_seq
        seq = Seq(whole_seq, generic_dna)
        seq_record = SeqRecord(seq)

        feature_1_loc = FeatureLocation(
                len(before), len(before) + len(feature_1_seq), strand=1)
        feature_1 = SeqFeature(feature_1_loc, type='CDS', id=1)
        seq_record.features.append(feature_1)

        # Simple table for testing.
        AA_TO_CODON_LIST_DICT = {
                'M': {
                    'ATG': {'usage': 1.0},
                },
                'T': {
                    'TTT': {'usage': 0.1},
                    'AAA': {'usage': 0.9},
                },
                'R': {
                    'CCC': {'usage': 0.8},
                    'GGG': {'usage': 0.1},
                    'TTT': {'usage': 0.1},
                }
        }
        CODON_USAGE_MEMEX = CodonUsageMemex(AA_TO_CODON_LIST_DICT)

        kwargs = {
                'original_codon_usage_memex': CODON_USAGE_MEMEX,
                'refactored_codon_usage_memex': CODON_USAGE_MEMEX
        }

        profile = CodonRarityFeatureProfile(
                feature_1.id, seq_record, **kwargs)

        EXPECTED_PROFILE_VALUES = [
                (0.9 + 0.8 + 0.1 + 0.1 + 1.0) / 5,
                (0.8 + 0.1 + 0.1 + 1.0 + 0.1) / 5,
                (0.1 + 0.1 + 1.0 + 0.1 + 0.1) / 5,
        ]
        self.assertEqual(EXPECTED_PROFILE_VALUES, profile.values)




if __name__ == '__main__':
    unittest.main()
